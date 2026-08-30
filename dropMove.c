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
#include <float.h>
#include <limits.h>
#include <sys/stat.h>

#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "src-local/activity.h"
#include "src-local/delay-history.h"

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
double delay_time = 0., delay_sample_interval = HUGE;
int min_level = 5, max_level = 9, drops_x = 1, drops_y = 1;
int delay_samples = 16;
int resume_run = 0, exit_status = 0, skip_initial_delay_properties = 0;
char output_dir[256] = "intermediate";
char restart_file[512] = "";

scalar * concentration_history = NULL;
DelayHistoryClock concentration_history_clock = {NULL, 0, 0, 0};

static void usage (const char * program)
{
  if (pid() == 0)
    fprintf (stderr,
      "Usage: %s [--params FILE] [key=value ...]\n"
      "Keys: Oh Ca Pe AcNum domain max_level min_level tmax tsnap\n"
      "      drops drops_x drops_y drop_radius drop_spacing ke_limit\n"
      "      movement_threshold t0 delay_samples output_dir\n"
      "      restart_file resume\n",
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

static int parse_delay_time (const char * value, double * parsed)
{
  char * end = NULL;
  errno = 0;
  double candidate = strtod (value, &end);
  while (end && isspace ((unsigned char) *end))
    end++;
  if (errno || end == value || !end || *end || !isfinite(candidate))
    return 0;
  *parsed = candidate;
  return 1;
}

static int parse_delay_samples (const char * value, int * parsed)
{
  char * end = NULL;
  errno = 0;
  long candidate = strtol (value, &end, 10);
  while (end && isspace ((unsigned char) *end))
    end++;
  if (errno || end == value || !end || *end ||
      candidate < INT_MIN || candidate > INT_MAX)
    return 0;
  *parsed = (int) candidate;
  return 1;
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
  else if (!strcmp (key, "t0")) {
    if (!parse_delay_time (value, &delay_time)) {
      if (pid() == 0)
        fprintf (stderr, "Invalid t0 value '%s'.\n", value);
      return 0;
    }
  }
  else if (!strcmp (key, "delay_samples")) {
    if (!parse_delay_samples (value, &delay_samples)) {
      if (pid() == 0)
        fprintf (stderr, "Invalid delay_samples value '%s'.\n", value);
      return 0;
    }
  }
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
    if (!strchr (line, '\n') && !feof(input)) {
      if (pid() == 0)
        fprintf (stderr, "Parameter line is too long at %s:%d.\n",
                 path, line_number);
      fclose (input);
      return 0;
    }
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
    int assignment_length =
      snprintf (assignment, sizeof(assignment), "%s", argv[arg]);
    if (assignment_length < 0 ||
        (size_t) assignment_length >= sizeof(assignment)) {
      if (pid() == 0)
        fprintf (stderr, "Argument is too long.\n");
      return 0;
    }
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
      !isfinite(delay_time) || delay_time < 0. ||
      delay_samples < 2 || delay_samples > 256 ||
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

static double delay_time_tolerance (double first, double second)
{
  return 64.*DBL_EPSILON*max(1., max(fabs(first), fabs(second)));
}

static int delay_metadata_path (const char * checkpoint, char * path,
                                size_t path_size)
{
  int length = snprintf (path, path_size, "%s.delay-history", checkpoint);
  return length > 0 && (size_t) length < path_size;
}

static double concentration_time_before_step (double event_time)
{
  return event_time > delay_time_tolerance (event_time, 0.) ?
    event_time - 0.5*dt : 0.;
}

/**
The default Basilisk restore path silently skips fields present only in the
dump. Inspect the name manifest first so an instantaneous run cannot consume a
delayed checkpoint and a delayed run cannot accept zero-filled history.
*/
static int checkpoint_delay_layout_matches (const char * checkpoint)
{
  int ok = 1;
  if (pid() == 0) {
    FILE * input = fopen (checkpoint, "r");
    int expected = delay_time > 0. ? delay_samples + 2 : 0;
    struct DumpHeader header = {0};

    if (!input || fread (&header, sizeof(header), 1, input) != 1 ||
        header.len < 1 || header.len > 100000) {
      fprintf (stderr, "Cannot inspect checkpoint field manifest %s.\n",
               checkpoint);
      ok = 0;
    }
    else if (header.version == 161020) {
      if (expected) {
        fprintf (stderr,
                 "Positive-delay restart requires a named-field checkpoint.\n");
        ok = 0;
      }
    }
    else if (header.version != dump_version) {
      fprintf (stderr, "Checkpoint version does not match this Basilisk.\n");
      ok = 0;
    }
    else
      ok = delay_history_manifest_matches (input, header.len, expected);

    if (input)
      fclose (input);
    if (!ok)
      fprintf (stderr,
               "Checkpoint delay fields do not match t0 and delay_samples.\n");
  }
#if _MPI
  MPI_Bcast (&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  return ok;
}

static int initialise_concentration_history (void)
{
  if (delay_time == 0.) {
    delay_sample_interval = HUGE;
    return 1;
  }

  delay_sample_interval = delay_time/delay_samples;
  if (!isfinite(delay_sample_interval) || delay_sample_interval <= 0.) {
    if (pid() == 0)
      fprintf (stderr, "t0/delay_samples is not a usable time interval.\n");
    return 0;
  }

  concentration_history_clock.capacity = delay_samples + 2;
  concentration_history_clock.times =
    (double *) malloc (concentration_history_clock.capacity*sizeof(double));
  if (!concentration_history_clock.times) {
    if (pid() == 0)
      fprintf (stderr, "Cannot allocate delay-history timestamps.\n");
    return 0;
  }
  delay_history_clock_reset (&concentration_history_clock);

  for (int slot_number = 0;
       slot_number < concentration_history_clock.capacity; slot_number++) {
    scalar slot = new scalar;
    char name[64];
    snprintf (name, sizeof(name), "cL_history_%04d", slot_number);
    free (slot.name);
    slot.name = strdup (name);
    scalar_clone (slot, cL);
#if TREE
    slot.refine = refine_bilinear;
    set_prolongation (slot, refine_bilinear);
    set_restriction (slot, restriction_volume_average);
#endif
    concentration_history = list_append (concentration_history, slot);
  }
  reset (concentration_history, 0.);
  return 1;
}

static int record_concentration_history (double sample_time)
{
  int newest = delay_history_newest (&concentration_history_clock);
  if (newest >= 0) {
    double latest = concentration_history_clock.times[newest];
    double tolerance = delay_time_tolerance (sample_time, latest);
    if (fabs(sample_time - latest) <= tolerance)
      return 1;
    if (sample_time < latest) {
      if (pid() == 0)
        fprintf (stderr,
                 "Delay-history sample time moved backwards: %g < %g.\n",
                 sample_time, latest);
      return 0;
    }
  }

  int slot_number = concentration_history_clock.next;
  scalar slot = concentration_history[slot_number];
  foreach()
    slot[] = cL[];
  boundary ({slot});
  return delay_history_record (&concentration_history_clock, sample_time) ==
    slot_number;
}

static int write_delay_metadata (const char * checkpoint)
{
  int ok = 1;
  if (pid() == 0) {
    char path[1024] = "", temporary[1032] = "";
    int temporary_length = -1;
    if (!delay_metadata_path (checkpoint, path, sizeof(path))) {
      fprintf (stderr, "Delay-history metadata path is too long.\n");
      ok = 0;
    }
    if (ok && delay_time == 0.) {
      if (remove (path) != 0 && errno != ENOENT) {
        fprintf (stderr, "Cannot remove stale %s: %s.\n", path,
                 strerror(errno));
        ok = 0;
      }
    }
    if (ok && delay_time > 0. &&
        ((temporary_length = snprintf (temporary, sizeof(temporary),
                                       "%s~", path)) <= 0 ||
         (size_t) temporary_length >= sizeof(temporary))) {
      fprintf (stderr, "Delay-history metadata path is too long.\n");
      ok = 0;
    }
    FILE * output = ok && delay_time > 0. ? fopen (temporary, "w") : NULL;
    if (ok && !output) {
      if (delay_time > 0.) {
        fprintf (stderr, "Cannot write %s: %s.\n", temporary,
                 strerror(errno));
        ok = 0;
      }
    }
    if (ok && delay_time > 0.) {
      fprintf (output, "version 2\n");
      fprintf (output, "checkpoint_time %.17g\n", t);
      fprintf (output, "checkpoint_concentration_time %.17g\n",
               concentration_time_before_step (t));
      fprintf (output, "t0 %.17g\n", delay_time);
      fprintf (output, "delay_samples %d\n", delay_samples);
      fprintf (output, "sample_interval %.17g\n", delay_sample_interval);
      fprintf (output, "capacity %d\n",
               concentration_history_clock.capacity);
      fprintf (output, "count %d\n", concentration_history_clock.count);
      fprintf (output, "next %d\n", concentration_history_clock.next);
      for (int slot = 0; slot < concentration_history_clock.capacity; slot++)
        fprintf (output, "time %d %.17g\n", slot,
                 concentration_history_clock.times[slot]);
      int write_failed = ferror(output);
      int close_failed = fclose(output) != 0;
      if (write_failed || close_failed) {
        fprintf (stderr, "Cannot finish writing %s: %s.\n", temporary,
                 strerror(errno));
        ok = 0;
      }
      else if (rename (temporary, path) != 0) {
        fprintf (stderr, "Cannot publish %s: %s.\n", path, strerror(errno));
        ok = 0;
      }
    }
    if (!ok && temporary[0])
      remove (temporary);
  }
#if _MPI
  MPI_Bcast (&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  return ok;
}

static int read_delay_metadata (const char * checkpoint)
{
  int ok = 1;
  if (pid() == 0) {
    char path[1024] = "", key[32];
    FILE * input = NULL;
    double checkpoint_time = NAN, checkpoint_concentration_time = NAN;
    double stored_delay = NAN;
    double stored_interval = NAN;
    int version = 0, stored_samples = 0, stored_capacity = 0;
    int stored_count = 0, stored_next = 0;
    double * stored_times = NULL;

    if (!delay_metadata_path (checkpoint, path, sizeof(path)) ||
        !(input = fopen (path, "r"))) {
      fprintf (stderr,
               "Positive-delay restart requires companion metadata %s.\n",
               path);
      ok = 0;
    }
    if (ok &&
        (fscanf (input, "%31s %d", key, &version) != 2 ||
         strcmp (key, "version") || version != 2 ||
         fscanf (input, "%31s %lf", key, &checkpoint_time) != 2 ||
         strcmp (key, "checkpoint_time") ||
         fscanf (input, "%31s %lf", key,
                 &checkpoint_concentration_time) != 2 ||
         strcmp (key, "checkpoint_concentration_time") ||
         fscanf (input, "%31s %lf", key, &stored_delay) != 2 ||
         strcmp (key, "t0") ||
         fscanf (input, "%31s %d", key, &stored_samples) != 2 ||
         strcmp (key, "delay_samples") ||
         fscanf (input, "%31s %lf", key, &stored_interval) != 2 ||
         strcmp (key, "sample_interval") ||
         fscanf (input, "%31s %d", key, &stored_capacity) != 2 ||
         strcmp (key, "capacity") ||
         fscanf (input, "%31s %d", key, &stored_count) != 2 ||
         strcmp (key, "count") ||
         fscanf (input, "%31s %d", key, &stored_next) != 2 ||
         strcmp (key, "next")))
      ok = 0;

    if (ok &&
        (!delay_history_metadata_values_valid (
           checkpoint_time, checkpoint_concentration_time, stored_delay,
           stored_interval,
           delay_time_tolerance (checkpoint_concentration_time,
                                 checkpoint_time)) ||
         stored_samples != delay_samples ||
         stored_capacity != concentration_history_clock.capacity ||
         stored_count < 1 || stored_count > stored_capacity ||
         stored_next < 0 || stored_next >= stored_capacity ||
         (stored_count < stored_capacity && stored_next != stored_count) ||
         fabs(stored_delay - delay_time) >
           delay_time_tolerance (stored_delay, delay_time) ||
         fabs(stored_interval - delay_sample_interval) >
           delay_time_tolerance (stored_interval, delay_sample_interval) ||
         fabs(checkpoint_time - t) >
           delay_time_tolerance (checkpoint_time, t)))
      ok = 0;

    if (ok) {
      stored_times = (double *) malloc (stored_capacity*sizeof(double));
      ok = stored_times != NULL;
    }
    for (int slot = 0; ok && slot < stored_capacity; slot++) {
      int stored_slot = -1;
      if (fscanf (input, "%31s %d %lf", key, &stored_slot,
                  &stored_times[slot]) != 3 || strcmp (key, "time") ||
          stored_slot != slot)
        ok = 0;
    }
    if (input)
      fclose (input);

    DelayHistoryClock stored_clock = {
      stored_times, stored_capacity, stored_count, stored_next
    };
    if (ok) {
      int oldest = delay_history_oldest (&stored_clock);
      double previous = stored_clock.times[oldest];
      if (!isfinite(previous))
        ok = 0;
      for (int offset = 1; ok && offset < stored_count; offset++) {
        int slot = (oldest + offset)%stored_capacity;
        double current = stored_clock.times[slot];
        if (!isfinite(current) || current <= previous)
          ok = 0;
        previous = current;
      }
      DelayHistoryBracket bracket;
      if (ok && !delay_history_bracket (
            &stored_clock, checkpoint_concentration_time - delay_time,
            delay_time_tolerance (checkpoint_concentration_time,
                                  delay_time), &bracket))
        ok = 0;
    }

    if (ok) {
      concentration_history_clock.count = stored_count;
      concentration_history_clock.next = stored_next;
      for (int slot = 0; slot < stored_capacity; slot++)
        concentration_history_clock.times[slot] = stored_times[slot];
    }
    else
      fprintf (stderr,
               "Delay-history metadata is invalid or does not match t0, "
               "delay_samples and checkpoint time.\n");
    free (stored_times);
  }
#if _MPI
  MPI_Bcast (&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ok) {
    MPI_Bcast (&concentration_history_clock.count, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    MPI_Bcast (&concentration_history_clock.next, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    MPI_Bcast (concentration_history_clock.times,
               concentration_history_clock.capacity, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
  }
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
  if (!initialise_concentration_history())
    return 2;
  skip_initial_delay_properties = delay_time > 0.;
  run();
  return exit_status;
}

event init (i = 0)
{
  if (resume_run && !checkpoint_delay_layout_matches (restart_file)) {
    exit_status = 2;
    return 1;
  }
  int restored = resume_run && restore (file = restart_file);
  if (restored) {
    if (delay_time > 0. && !read_delay_metadata (restart_file)) {
      exit_status = 2;
      return 1;
    }
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
    if (delay_time > 0.) {
      delay_history_clock_reset (&concentration_history_clock);
      if (!record_concentration_history (0.)) {
        exit_status = 2;
        return 1;
      }
    }
  }
}

event delay_history_sample (t += delay_sample_interval; t <= end_time)
{
  double sample_time = concentration_time_before_step (t);
  if (delay_time > 0. && !record_concentration_history (sample_time)) {
    exit_status = 2;
    return 1;
  }
}

event properties (i++)
{
  if (delay_time == 0.) {
    foreach()
      sigmaf[] = 1./ca + 4.*cL[];
  }
  else {
    if (skip_initial_delay_properties) {
      skip_initial_delay_properties = 0;
      return 0;
    }
    DelayHistoryBracket bracket;
    double action_time = t + 0.5*dt;
    double target = action_time - delay_time;
    if (!delay_history_bracket (
          &concentration_history_clock, target,
          delay_time_tolerance (action_time, delay_time), &bracket)) {
      if (pid() == 0)
        fprintf (stderr,
                 "Delay history does not bracket tracer time - t0 = %.17g "
                 "at tracer time %.17g.\n", target, action_time);
      exit_status = 2;
      return 1;
    }
    scalar lower = concentration_history[bracket.lower];
    scalar upper = concentration_history[bracket.upper];
    foreach() {
      double delayed = delay_history_value (
        cL[], lower[], upper[], bracket.weight, delay_time);
      sigmaf[] = 1./ca + 4.*delayed;
    }
  }
}

event adapt (i++)
{
  foreach()
    KAPPA[] = distance_curvature (point, d);
  if (delay_time == 0.)
    adapt_wavelet ({f, u.x, u.y, cL, KAPPA},
                   (double[]){fraction_tolerance, velocity_tolerance,
                              velocity_tolerance, concentration_tolerance,
                              curvature_tolerance},
                   max_level, min_level);
  else {
    scalar * monitored = list_copy ({f, u.x, u.y, cL, KAPPA});
    int oldest = delay_history_oldest (&concentration_history_clock);
    for (int offset = 0; offset < concentration_history_clock.count;
         offset++) {
      int slot = (oldest + offset)%concentration_history_clock.capacity;
      monitored = list_append (monitored, concentration_history[slot]);
    }
    int fields = list_len (monitored);
    double * tolerances = (double *) malloc (fields*sizeof(double));
    if (!tolerances) {
      free (monitored);
      if (pid() == 0)
        fprintf (stderr, "Cannot allocate delay-history AMR tolerances.\n");
      exit_status = 2;
      return 1;
    }
    tolerances[0] = fraction_tolerance;
    tolerances[1] = velocity_tolerance;
    tolerances[2] = velocity_tolerance;
    tolerances[3] = concentration_tolerance;
    tolerances[4] = curvature_tolerance;
    for (int field = 5; field < fields; field++)
      tolerances[field] = concentration_tolerance;
    adapt_wavelet (monitored, tolerances, max_level, min_level);
    free (tolerances);
    free (monitored);
  }
}

event outputs (t = 0.; t += snapshot_interval; t <= end_time)
{
  char snapshot[512];
  snprintf (snapshot, sizeof(snapshot), "%s/snapshot-%012.6f",
            output_dir, t);
  dump (file = snapshot);
  if (!write_delay_metadata (snapshot)) {
    exit_status = 2;
    return 1;
  }
  dump (file = restart_file);
  if (!write_delay_metadata (restart_file)) {
    exit_status = 2;
    return 1;
  }
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
    if (!write_delay_metadata (failure_dump)) {
      exit_status = 2;
      return 1;
    }
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

event cleanup_delay_history (t = end)
{
  if (concentration_history) {
    delete (concentration_history);
    free (concentration_history);
    concentration_history = NULL;
  }
  free (concentration_history_clock.times);
  concentration_history_clock.times = NULL;
}
