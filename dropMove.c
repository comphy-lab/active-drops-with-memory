/**
 * Planar chemically fuelled active-drop simulation.
 *
 * The interface is represented with coupled level-set/VOF (CLSVOF);
 * concentration-dependent surface tension drives Marangoni flow.
 * Runtime parameters can be supplied as key=value arguments or through
 * a text file passed with --params FILE.
 */

#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>

#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "src-local/activity.h"

scalar cL[], *stracers = {cL};
scalar sigmaf[], KAPPA[];

cL[top] = dirichlet(0.);
cL[right] = dirichlet(0.);
cL[left] = dirichlet(0.);
cL[bottom] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

double oh = 1., ca = 0.1, pe = 1.6, activity = 1.;
double domain_size = 10., drop_radius = 1., drop_spacing = 3.;
double end_time = 50., snapshot_interval = 0.1;
double velocity_tolerance = 1e-3, fraction_tolerance = 1e-3;
double concentration_tolerance = 1e-3, curvature_tolerance = 1e-3;
double kinetic_energy_limit = 1e3;
double movement_threshold = -1.;
int min_level = 5, max_level = 9, drops_x = 1, drops_y = 1;
int resume_run = 0, exit_status = 0;
char output_dir[256] = "intermediate";
char restart_file[512] = "";

static void usage (const char * program)
{
  if (pid() == 0)
    fprintf (stderr,
      "Usage: %s [--params FILE] [key=value ...]\n"
      "Keys: Oh Ca Pe AcNum domain max_level min_level tmax tsnap\n"
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
  if (!strcmp (key, "Oh"))
    oh = atof (value);
  else if (!strcmp (key, "Ca"))
    ca = atof (value);
  else if (!strcmp (key, "Pe"))
    pe = atof (value);
  else if (!strcmp (key, "AcNum"))
    activity = atof (value);
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
    snprintf (assignment, sizeof(assignment), "%s", argv[arg]);
    if (strchr (assignment, '=')) {
      if (!parse_assignment (assignment))
        return 0;
    }
    else if (!legacy_pe_used) {
      pe = atof (assignment);
      legacy_pe_used = 1;
    }
    else {
      if (pid() == 0)
        fprintf (stderr, "Invalid argument '%s'.\n", argv[arg]);
      return 0;
    }
  }

  if (oh <= 0. || ca <= 0. || pe <= 0. || domain_size <= 0. ||
      drop_radius <= 0. || drop_spacing <= 0. || end_time <= 0. ||
      snapshot_interval <= 0. || kinetic_energy_limit <= 0. ||
      min_level < 1 || max_level < min_level || max_level > 20 ||
      drops_x < 1 || drops_y < 1 ||
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
  init_grid (1 << min_level);

  d.sigmaf = sigmaf;
  rho1 = 4./sq(oh);
  rho2 = 4./sq(oh);
  mu1 = 1.;
  mu2 = 1.;
  cL.inverse = true;
  cL.A = activity;
  cL.D = 1./pe;
  run();
  return exit_status;
}

event init (i = 0)
{
  int restored = resume_run && restore (file = restart_file);
  if (restored) {
    if (pid() == 0)
      fprintf (stderr, "Restarted from %s at t=%g, i=%d.\n",
               restart_file, t, i);
  }
  else if (resume_run) {
    if (pid() == 0)
      fprintf (stderr, "Cannot restore checkpoint %s.\n", restart_file);
    exit_status = 2;
    return 1;
  }
  else {
    refine (fabs(initial_distance (x, y)) < 2.*drop_radius &&
            level < max_level);
    foreach() {
      d[] = initial_distance (x, y);
      foreach_dimension()
        u.x[] = 0.;
      cL[] = 0.;
      sigmaf[] = 1./ca + 4.*cL[];
    }
  }
}

event properties (i++)
{
  foreach()
    sigmaf[] = 1./ca + 4.*cL[];
}

event adapt (i++)
{
  foreach()
    KAPPA[] = distance_curvature (point, d);
  adapt_wavelet ({f, u.x, u.y, cL, KAPPA},
                 (double[]){fraction_tolerance, velocity_tolerance,
                            velocity_tolerance, concentration_tolerance,
                            curvature_tolerance},
                 max_level, min_level);
}

event outputs (t = 0.; t += snapshot_interval; t <= end_time)
{
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

  if (pid() == 0) {
    static FILE * log = NULL;
    char log_path[512];
    snprintf (log_path, sizeof(log_path), "%s/log.dat", output_dir);
    if (!log) {
      log = fopen (log_path, resume_run ? "a" : "w");
      if (log && !resume_run)
        fprintf (log, "i t ke\n");
    }
    fprintf (stderr, "%d %g %.8e\n", i, t, ke);
    if (log) {
      fprintf (log, "%d %g %.8e\n", i, t, ke);
      fflush (log);
    }
  }

  if (!isfinite(ke) || ke >= kinetic_energy_limit) {
    char failure_dump[512];
    snprintf (failure_dump, sizeof(failure_dump), "%s/failure-%d",
              output_dir, i);
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
