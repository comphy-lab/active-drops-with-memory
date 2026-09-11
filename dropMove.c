/**
# Deformable drops consuming a finite exterior fuel reservoir

Full incompressible Navier--Stokes flow couples to phase-confined fuel
transport and saturating interfacial uptake. Depletion raises surface
tension. Both spatial directions are periodic. See the README for the
dimensionless model, numerical approximations and biological limitations.
*/

#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>

#include "navier-stokes/centered.h"
static bool fuel_repair_roundoff (void);
extern int exit_status;
/* Register before vof.h: the inherited core sweep runs first. */
event vof (i++) {
  if (!fuel_repair_roundoff()) {
    fprintf(stderr, "Fuel VOF residue exceeds the conservative roundoff budget.\n");
    exit_status = 3;
    return 1;
  }
}
#include "vof.h"
static bool initialise_fuel (void);
/* Inherited events run in reverse order: user distance, CLSVOF fraction,
   fuel initialisation, then centered fluid initialisation. */
event init (i = 0) {
  if (!initialise_fuel()) return 1;
}
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "src-local/finite-fuel.h"

scalar sigmaf[], KAPPA[];
double reynolds = 0.01, ca = 0.1, pe = 1.6, damkohler = 1.;
double gamma_slope = 4., initial_fuel = 1., half_saturation = 0.1;
double viscosity_ratio = 1., density_ratio = 1.;
double domain_size = 10., drop_radius = 1., drop_spacing = 3.;
double end_time = 50., snapshot_interval = 0.1;
double velocity_tolerance = 1e-3, fraction_tolerance = 1e-3;
double concentration_tolerance = 1e-3, curvature_tolerance = 1e-3;
double kinetic_energy_limit = 1e3;
double movement_threshold = -1.;
int min_level = 5, max_level = 9, drops_x = 1, drops_y = 1;
int resume_run = 0, exit_status = 0;
static bool resume_pending = true;
char output_dir[256] = "intermediate";
char restart_file[512] = "";

static void usage (const char * program)
{
  if (pid() == 0)
    fprintf (stderr,
      "Usage: %s [--params FILE] [key=value ...]\n"
      "Keys: Re Ca Pe Da GammaSlope InitialFuel FuelHalfSaturation\n"
      "      viscosityRatio densityRatio domain max_level min_level tmax tsnap\n"
      "      drops drops_x drops_y drop_radius drop_spacing ke_limit\n"
      "      movement_threshold output_dir restart_file resume\n",
      program);
}

static char * trim (char * text)
{
  while (isspace ((unsigned char) *text))
    text++;
  char * end = text + strlen (text);
  while (end > text && isspace ((unsigned char) end[-1]))
    *--end = '\0';
  return text;
}

static int set_parameter (const char * key, const char * value)
{
  if (!strcmp(key, "Oh") || !strcmp(key, "AcNum")) {
    fprintf(stderr, "Legacy %s is not a finite-fuel parameter: use Re, Da and GammaSlope; Oh=sqrt(Ca/Re).\n", key);
    return 0;
  }
  if (!strcmp(key, "output_dir") || !strcmp(key, "restart_file")) {
    char * destination = !strcmp(key, "output_dir") ? output_dir : restart_file;
    size_t capacity = !strcmp(key, "output_dir") ? sizeof(output_dir) : sizeof(restart_file);
    if (!*value || strlen(value) >= capacity) return 0;
    snprintf(destination, capacity, "%s", value);
    return 1;
  }
  char * end;
  errno = 0;
  double number = strtod(value, &end);
  if (errno || end == value || *end || !isfinite(number)) return 0;
  bool integer = !strcmp(key, "max_level") || !strcmp(key, "min_level") ||
    !strcmp(key, "drops") || !strcmp(key, "drops_x") || !strcmp(key, "drops_y") ||
    !strcmp(key, "resume");
  if (integer && (number != floor(number) || number < 0. || number > 1000000.)) return 0;
  if (!strcmp (key, "Re"))
    reynolds = number;
  else if (!strcmp (key, "Ca"))
    ca = number;
  else if (!strcmp (key, "Pe"))
    pe = number;
  else if (!strcmp(key, "Da")) damkohler = number;
  else if (!strcmp(key, "GammaSlope")) gamma_slope = number;
  else if (!strcmp(key, "InitialFuel")) initial_fuel = number;
  else if (!strcmp(key, "FuelHalfSaturation")) half_saturation = number;
  else if (!strcmp(key, "viscosityRatio")) viscosity_ratio = number;
  else if (!strcmp(key, "densityRatio")) density_ratio = number;
  else if (!strcmp (key, "domain"))
    domain_size = atof (value);
  else if (!strcmp (key, "max_level"))
    max_level = atoi (value);
  else if (!strcmp (key, "min_level"))
    min_level = atoi (value);
  else if (!strcmp (key, "tmax"))
    end_time = atof (value);
  else if (!strcmp (key, "tsnap"))
    snapshot_interval = atof (value);
  else if (!strcmp (key, "drops")) {
    int drops = atoi (value);
    if (drops == 1)
      drops_x = drops_y = 1;
    else if (drops == 2)
      drops_x = 2, drops_y = 1;
    else if (drops == 4)
      drops_x = drops_y = 2;
    else
      return 0;
  }
  else if (!strcmp (key, "drops_x"))
    drops_x = atoi (value);
  else if (!strcmp (key, "drops_y"))
    drops_y = atoi (value);
  else if (!strcmp (key, "drop_radius"))
    drop_radius = atof (value);
  else if (!strcmp (key, "drop_spacing"))
    drop_spacing = atof (value);
  else if (!strcmp (key, "ke_limit"))
    kinetic_energy_limit = atof (value);
  else if (!strcmp (key, "movement_threshold"))
    movement_threshold = atof (value);
  else if (!strcmp (key, "output_dir"))
    snprintf (output_dir, sizeof(output_dir), "%s", value);
  else if (!strcmp (key, "restart_file"))
    snprintf (restart_file, sizeof(restart_file), "%s", value);
  else if (!strcmp (key, "resume"))
    resume_run = atoi (value);
  else {
    if (pid() == 0)
      fprintf (stderr, "Unknown parameter '%s'.\n", key);
    return 0;
  }
  return 1;
}

static int parse_assignment (char * assignment)
{
  char * equals = strchr (assignment, '=');
  if (!equals)
    return 0;
  *equals = '\0';
  return set_parameter (trim (assignment), trim (equals + 1));
}

static int read_parameter_file (const char * path)
{
  FILE * input = fopen (path, "r");
  if (!input) {
    if (pid() == 0)
      fprintf (stderr, "Cannot open parameter file '%s': %s\n",
               path, strerror (errno));
    return 0;
  }
  char line[1024];
  int line_number = 0;
  while (fgets (line, sizeof(line), input)) {
    line_number++;
    char * comment = strchr (line, '#');
    if (comment)
      *comment = '\0';
    char * content = trim (line);
    if (!*content)
      continue;
    if (!parse_assignment (content)) {
      if (pid() == 0)
        fprintf (stderr, "Invalid parameter at %s:%d.\n", path, line_number);
      fclose (input);
      return 0;
    }
  }
  fclose (input);
  return 1;
}

static int parse_parameters (int argc, char ** argv)
{
  int legacy_pe_used = 0;
  for (int arg = 1; arg < argc; arg++) {
    if (!strcmp (argv[arg], "--help")) {
      usage (argv[0]);
      return 0;
    }
    if (!strcmp (argv[arg], "--params")) {
      if (++arg == argc || !read_parameter_file (argv[arg]))
        return 0;
      continue;
    }
    if (!strncmp (argv[arg], "--params=", 9)) {
      if (!read_parameter_file (argv[arg] + 9))
        return 0;
      continue;
    }
    char assignment[1024];
    if (strlen(argv[arg]) >= sizeof(assignment)) return 0;
    snprintf (assignment, sizeof(assignment), "%s", argv[arg]);
    if (strchr (assignment, '=')) {
      if (!parse_assignment (assignment))
        return 0;
    }
    else if (!legacy_pe_used) {
      if (!set_parameter("Pe", assignment)) return 0;
      legacy_pe_used = 1;
    }
    else {
      if (pid() == 0)
        fprintf (stderr, "Invalid argument '%s'.\n", argv[arg]);
      return 0;
    }
  }

  if (reynolds <= 0. || ca <= 0. || pe <= 0. || domain_size <= 0. ||
      damkohler < 0. || gamma_slope < 0. || initial_fuel < 0. || initial_fuel > 1. ||
      half_saturation <= 0. || viscosity_ratio <= 0. || density_ratio <= 0. ||
      !isfinite(reynolds*density_ratio) || !isfinite(1./ca + gamma_slope) ||
      !isfinite(1./pe) || !isfinite(damkohler/pe) ||
      fraction_tolerance <= 0. || velocity_tolerance <= 0. ||
      concentration_tolerance <= 0. || curvature_tolerance <= 0. ||
      (resume_run != 0 && resume_run != 1) ||
      drop_radius <= 0. || drop_spacing <= 0. || end_time <= 0. ||
      snapshot_interval <= 0. || kinetic_energy_limit <= 0. ||
      min_level < 1 || max_level < min_level || max_level > 20 ||
      drops_x < 1 || drops_y < 1 ||
      ((drops_x > 1 || drops_y > 1) && drop_spacing <= 2.*drop_radius) ||
      (drops_x - 1)*drop_spacing + 2.*drop_radius >= domain_size ||
      (drops_y - 1)*drop_spacing + 2.*drop_radius >= domain_size) {
    if (pid() == 0)
      fprintf (stderr, "Invalid parameter values or drop lattice.\n");
    return 0;
  }
  return 1;
}

static int make_directory_tree (const char * path)
{
  char directory[256];
  snprintf (directory, sizeof(directory), "%s", path);
  for (char * separator = directory + 1; *separator; separator++)
    if (*separator == '/') {
      *separator = '\0';
      if (mkdir (directory, 0775) != 0 && errno != EEXIST)
        return 0;
      *separator = '/';
    }
  if (mkdir (directory, 0775) == 0)
    return 1;
  if (errno != EEXIST)
    return 0;
  struct stat status;
  return stat (directory, &status) == 0 && S_ISDIR (status.st_mode);
}

static int prepare_output_directory (void)
{
  int ok = 1;
  if (pid() == 0 && !make_directory_tree (output_dir)) {
    fprintf (stderr, "Cannot create output directory '%s': %s\n",
             output_dir, strerror (errno));
    ok = 0;
  }
#if _MPI
  MPI_Bcast (&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  return ok;
}

static double initial_distance (double px, double py)
{
  double signed_distance = -HUGE;
  for (int ix = 0; ix < drops_x; ix++)
    for (int iy = 0; iy < drops_y; iy++) {
      double centre_x = (ix - 0.5*(drops_x - 1))*drop_spacing;
      double centre_y = (iy - 0.5*(drops_y - 1))*drop_spacing;
      signed_distance =
        max (signed_distance,
             drop_radius -
             sqrt (sq(px - centre_x) + sq(py - centre_y)));
    }
  return signed_distance;
}

/**
Checkpoint metadata is stored inside the dump, so fuel accounting and
physics cannot become detached from the state. Legacy dumps lack the
version marker. Physical parameters and initial layout must match exactly.
*/
scalar meta_version[], meta_re[], meta_ca[], meta_pe[], meta_da[], meta_gamma[];
scalar meta_half[], meta_initial[], meta_mu[], meta_rho[], meta_domain[];
scalar meta_radius[], meta_spacing[], meta_nx[], meta_ny[];
scalar meta_inventory[], meta_consumed[], meta_roundoff[];
scalar meta_dt[], checkpoint_ux_low[], checkpoint_ux_high[];
scalar checkpoint_uy_low[], checkpoint_uy_high[];

static bool fuel_metadata (bool write)
{
  double values[] = {26091103., reynolds, ca, pe, damkohler, gamma_slope,
    half_saturation, initial_fuel, viscosity_ratio, density_ratio, domain_size,
    drop_radius, drop_spacing, drops_x, drops_y};
  scalar * fields = {meta_version, meta_re, meta_ca, meta_pe, meta_da, meta_gamma,
    meta_half, meta_initial, meta_mu, meta_rho, meta_domain,
    meta_radius, meta_spacing, meta_nx, meta_ny};
  int index = 0;
  for (scalar field in fields) {
    double expected = values[index++];
    if (write) {
      foreach() field[] = expected;
    }
    else {
      stats range = statsf(field);
      if (range.min != expected || range.max != expected) {
        fprintf(stderr, "Missing or incompatible finite-fuel checkpoint metadata: %s.\n", field.name);
        return false;
      }
    }
  }
  if (write) {
    foreach() {
      meta_inventory[] = fuel_initial_inventory;
      meta_consumed[] = fuel_consumed;
      meta_roundoff[] = fuel_roundoff_inventory;
      meta_dt[] = dt;
      checkpoint_ux_low[] = uf.x[];
      checkpoint_ux_high[] = uf.x[1];
      checkpoint_uy_low[] = uf.y[];
      checkpoint_uy_high[] = uf.y[0,1];
    }
  }
  else {
    stats initial = statsf(meta_inventory), consumed = statsf(meta_consumed);
    if (!isfinite(initial.min) || initial.min < 0. || initial.min != initial.max ||
        !isfinite(consumed.min) || consumed.min < 0. || consumed.min != consumed.max)
      return false;
    fuel_initial_inventory = initial.min;
    fuel_consumed = consumed.min;
    stats roundoff = statsf(meta_roundoff);
    if (!isfinite(roundoff.min) || roundoff.min < 0. || roundoff.min != roundoff.max ||
        roundoff.min > 4096.*DBL_EPSILON*fuel_initial_inventory) return false;
    fuel_roundoff_inventory = roundoff.min;
    stats timestep = statsf(meta_dt);
    if (!isfinite(timestep.min) || timestep.min <= 0. || timestep.min != timestep.max)
      return false;
    dt = timestep.min;
  }
  return true;
}

/**
Outputs save the state at the beginning of a physical timestep. Basilisk's
restore catch-up advances recurring iteration events past that iteration.
Re-arm those events so the normal solver advances the saved state once;
time-based outputs remain advanced and are not replayed.
*/
static void resume_iteration_events()
{
  for (Event * ev = Events; !ev->last; ev++)
    if (ev->expr[2] && ev->i == iter + 1 && ev->t < 0.)
      ev->i = iter;
}

static bool initialise_fuel (void)
{
  foreach()
    fuel[] = (1. - f[])*initial_fuel;
  fuel_consumed = 0.;
  fuel_roundoff_inventory = 0.;
  fuel_initial_inventory = fuel_inventory();
  p.nodump = pf.nodump = false;
  fuel_ready = true;
  if (!fuel_recover()) { exit_status = 3; return false; }
  return true;
}

int main (int argc, char ** argv)
{
  if (!parse_parameters (argc, argv))
    return 2;
  if (!prepare_output_directory())
    return 2;
  if (!restart_file[0])
    snprintf (restart_file, sizeof(restart_file), "%s/restart", output_dir);

  size (domain_size);
  origin (-0.5*domain_size, -0.5*domain_size);
  periodic(right);
  periodic(top);
  init_grid (1 << min_level);

  d.sigmaf = sigmaf;
  rho1 = reynolds*density_ratio;
  rho2 = reynolds;
  mu1 = viscosity_ratio;
  mu2 = 1.;
  fuel.inverse = true;
  fuel.gradient = zero;
  f.tracers = list_append(f.tracers, fuel);
  fuel_pe = pe; fuel_da = damkohler; fuel_half_saturation = half_saturation;
  fuel_roundoff_scale = initial_fuel;
  TOLERANCE = 1e-10;
  fprintf(stderr, "Re=%g Ca=%g Pe=%g Da=%g GammaSlope=%g Oh_derived=%g\n",
          reynolds, ca, pe, damkohler, gamma_slope, sqrt(ca/reynolds));
  run();
  return exit_status;
}

event init (i = 0)
{
    refine (fabs(initial_distance (x, y)) < 2.*drop_radius &&
            level < max_level);
    foreach() {
      d[] = initial_distance (x, y);
      foreach_dimension()
        u.x[] = 0.;
      fuel[] = 0.;
      sigmaf[] = 1./ca + gamma_slope*(1. - initial_fuel);
    }
}

/**
Restore after every inherited init action, so neither CLSVOF nor centered
initialisation rewrites the saved beginning-of-step state. Face fields
are represented by both sides of each leaf in the checkpoint; choosing
the fine leaf at a coarse/fine boundary retains its subface flux.
*/
event resume_state (i++) {
  if (!resume_pending || !resume_run) return 0;
  resume_pending = false;
  if (!restore(file = restart_file) || !fuel_metadata(false) || !fuel_recover() ||
      fabs(fuel_inventory() + fuel_consumed - fuel_initial_inventory) >
      1e-8*max(1., fuel_initial_inventory)) {
    fprintf(stderr, "Cannot resume a legacy, incompatible or unbalanced fuel checkpoint.\n");
    exit_status = 2;
    return 1;
  }
  foreach_face(x)
    uf.x[] = is_leaf(cell) ? checkpoint_ux_low[] : checkpoint_ux_high[-1];
  foreach_face(y)
    uf.y[] = is_leaf(cell) ? checkpoint_uy_low[] : checkpoint_uy_high[0,-1];
  boundary({uf});
  /* Refresh material face coefficients from the saved filtered fraction,
     without invoking the CLSVOF redistancing properties event. */
  foreach_face() {
    double fraction = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(fraction);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(fraction);
  }
  foreach()
    rhov[] = cm[]*rho(sf[]);
  resume_iteration_events();
}

event properties (i++)
{
  if (!fuel_ready) return 0;
  if (!fuel_recover()) { exit_status = 3; return 1; }
  fuel_surface_extension();
  int invalid = 0;
  foreach(reduction(max:invalid)) {
    sigmaf[] = 1./ca + gamma_slope*(1. - surface_fuel[]);
    if (!isfinite(sigmaf[]) || sigmaf[] <= 0.) invalid = 1;
  }
  if (invalid) {
    fprintf(stderr, "Nonpositive or nonfinite surface tension.\n");
    exit_status = 3;
    return 1;
  }
}

event tracer_diffusion (i++) {
  if (fuel_ready && (!fuel_diffuse(dt) || !fuel_consume(dt))) {
    fprintf(stderr, "Fuel positivity or inventory contract failed.\n");
    exit_status = 3;
    return 1;
  }
}

event adapt (i++)
{
  foreach()
    KAPPA[] = distance_curvature (point, d);
  double before = fuel_inventory();
  adapt_wavelet ({f, u.x, u.y, fuel, KAPPA},
                 (double[]){fraction_tolerance, velocity_tolerance,
                            velocity_tolerance, concentration_tolerance,
                            curvature_tolerance},
                 max_level, min_level);
  if (!fuel_recover() || fabs(fuel_inventory() - before) > 1e-8*max(1., before)) {
    fprintf(stderr, "Fuel adaptation contract failed.\n");
    exit_status = 3;
    return 1;
  }
}

event outputs (t = 0.; t += snapshot_interval; t <= end_time)
{
  fuel_metadata(true);
  char snapshot[512];
  snprintf (snapshot, sizeof(snapshot), "%s/snapshot-%012.6f",
            output_dir, t);
  dump (file = snapshot);
  dump (file = restart_file);
}

event logWriting (i++)
{
  double ke = 0.;
  double drop_area = 0., x_moment = 0., y_moment = 0.;
  foreach (reduction(+:ke) reduction(+:drop_area)
           reduction(+:x_moment) reduction(+:y_moment)) {
    ke += 0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
    double cell_area = clamp(f[], 0., 1.)*sq(Delta);
    drop_area += cell_area;
    x_moment += cell_area*x;
    y_moment += cell_area*y;
  }
  double displacement = drop_area > 0. ?
    sqrt (sq(x_moment/drop_area) + sq(y_moment/drop_area)) : 0.;
  double inventory = fuel_inventory();
  double balance = inventory + fuel_consumed - fuel_initial_inventory;

  if (pid() == 0) {
    static FILE * log = NULL;
    char log_path[512];
    snprintf (log_path, sizeof(log_path), "%s/log.dat", output_dir);
    if (!log) {
      log = fopen (log_path, resume_run ? "a" : "w");
      if (log && !resume_run)
        fprintf (log, "i t ke fuel_inventory fuel_consumed fuel_balance drop_area fuel_roundoff_transferred\n");
    }
    fprintf (stderr, "%d %g %.8e\n", i, t, ke);
    if (log) {
      fprintf (log, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
               i, t, ke, inventory, fuel_consumed, balance, drop_area, fuel_roundoff_inventory);
      fflush (log);
    }
  }

  if (!isfinite(ke) || ke >= kinetic_energy_limit || !isfinite(balance) ||
      fabs(balance) > 1e-7*max(1., fuel_initial_inventory) || !fuel_recover()) {
    char failure_dump[512];
    snprintf (failure_dump, sizeof(failure_dump), "%s/failure-%d",
              output_dir, i);
    fuel_metadata(true);
    dump (file = failure_dump);
    if (pid() == 0)
      fprintf (stderr,
               "Stopping cleanly: kinetic energy %.8e at i=%d, t=%g; "
               "state saved to %s.\n", ke, i, t, failure_dump);
    exit_status = 3;
    return 1;
  }
  if (movement_threshold > 0. && displacement >= movement_threshold) {
    if (pid() == 0)
      fprintf (stdout, "STATUS MOVED\n");
    return 1;
  }
}

event end (t = end_time)
{
  if (movement_threshold > 0. && pid() == 0)
    fprintf (stdout, "STATUS NOT_MOVED\n");
}
