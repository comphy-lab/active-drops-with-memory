#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../src-local/delay-history.h"

#define CHECK(condition) do {                                             \
  if (!(condition)) {                                                     \
    fprintf (stderr, "%s:%d: check failed: %s\n",                         \
             __FILE__, __LINE__, #condition);                             \
    exit (EXIT_FAILURE);                                                  \
  }                                                                       \
} while (0)

static void test_zero_delay_is_exact (void)
{
  double current = nextafter (1., 2.);
  CHECK (delay_history_value (current, -3., 9., 0.25, 0.) == current);
  CHECK (signbit (delay_history_value (-0., 1., 2., 0.5, 0.)));
}

static void test_startup_and_interpolation (void)
{
  double times[5];
  DelayHistoryClock clock = {times, 5, 0, 0};
  delay_history_clock_reset (&clock);
  CHECK (delay_history_record (&clock, 0.) == 0);
  CHECK (delay_history_record (&clock, 0.25) == 1);
  CHECK (delay_history_record (&clock, 0.5) == 2);

  DelayHistoryBracket bracket;
  CHECK (delay_history_bracket (&clock, -0.5, 1e-12, &bracket));
  CHECK (bracket.lower == 0 && bracket.upper == 0);

  CHECK (delay_history_bracket (&clock, 0.375, 1e-12, &bracket));
  CHECK (bracket.lower == 1 && bracket.upper == 2);
  CHECK (fabs(bracket.weight - 0.5) < 1e-14);
  CHECK (delay_history_value (99., 2., 6., bracket.weight, 1.) == 4.);

  CHECK (delay_history_bracket (&clock, 0.5, 1e-12, &bracket));
  CHECK (bracket.lower == 2 && bracket.upper == 2);
  CHECK (!delay_history_bracket (&clock, 0.6, 1e-12, &bracket));
}

static void test_wrapped_ring (void)
{
  double times[4];
  DelayHistoryClock clock = {times, 4, 0, 0};
  delay_history_clock_reset (&clock);
  for (int sample = 0; sample <= 4; sample++)
    CHECK (delay_history_record (&clock, sample) == sample%4);

  CHECK (clock.count == 4 && clock.next == 1);
  CHECK (delay_history_oldest (&clock) == 1);
  CHECK (delay_history_newest (&clock) == 0);

  DelayHistoryBracket bracket;
  CHECK (!delay_history_bracket (&clock, 0.5, 1e-12, &bracket));
  CHECK (delay_history_bracket (&clock, 2.5, 1e-12, &bracket));
  CHECK (bracket.lower == 2 && bracket.upper == 3);
  CHECK (fabs(bracket.weight - 0.5) < 1e-14);
}

static void write_manifest_name (FILE * output, const char * name)
{
  unsigned length = (unsigned) strlen (name);
  CHECK (fwrite (&length, sizeof(length), 1, output) == 1);
  CHECK (fwrite (name, sizeof(char), length, output) == length);
}

static FILE * manifest_stream (const char ** names, int count)
{
  FILE * output = tmpfile();
  CHECK (output);
  for (int index = 0; index < count; index++)
    write_manifest_name (output, names[index]);
  rewind (output);
  return output;
}

static void test_checkpoint_manifests (void)
{
  const char * no_history[] = {"size", "cL", "sigmaf"};
  FILE * input = manifest_stream (no_history, 3);
  CHECK (delay_history_manifest_matches (input, 3, 0));
  fclose (input);

  const char * complete[] = {
    "size", "cL", "cL_history_0000", "cL_history_0001",
    "cL_history_0002", "sigmaf"
  };
  input = manifest_stream (complete, 6);
  CHECK (delay_history_manifest_matches (input, 6, 3));
  fclose (input);

  input = manifest_stream (complete, 6);
  CHECK (!delay_history_manifest_matches (input, 6, 0));
  fclose (input);

  const char * missing[] = {
    "size", "cL", "cL_history_0000", "cL_history_0002", "sigmaf"
  };
  input = manifest_stream (missing, 5);
  CHECK (!delay_history_manifest_matches (input, 5, 3));
  fclose (input);

  const char * duplicate[] = {
    "size", "cL_history_0000", "cL_history_0000",
    "cL_history_0001", "cL_history_0002"
  };
  input = manifest_stream (duplicate, 5);
  CHECK (!delay_history_manifest_matches (input, 5, 3));
  fclose (input);
}

static void test_variable_timestep_stagger (void)
{
  enum { samples = 4, capacity = samples + 2 };
  double times[capacity];
  DelayHistoryClock clock = {times, capacity, 0, 0};
  delay_history_clock_reset (&clock);
  CHECK (delay_history_record (&clock, 0.) == 0);
  const double interval = 0.25, delay = samples*interval;
  const double previous_dt[] = {
    0.20, 0.08, 0.24, 0.12, 0.18, 0.06, 0.22, 0.10,
    0.16, 0.04, 0.20, 0.14, 0.08
  };
  const double current_dt[] = {
    0.08, 0.24, 0.12, 0.18, 0.06, 0.22, 0.10, 0.16,
    0.04, 0.20, 0.14, 0.08, 0.24
  };
  for (int step = 1; step <= 12; step++) {
    double nominal = step*interval;
    double stored_time = nominal - 0.5*previous_dt[step];
    CHECK (delay_history_record (&clock, stored_time) >= 0);
    double action_time = nominal + 0.5*current_dt[step];
    DelayHistoryBracket bracket;
    CHECK (delay_history_bracket (&clock, action_time - delay, 1e-12,
                                   &bracket));
  }
}

static void test_metadata_identity_values (void)
{
  CHECK (delay_history_metadata_values_valid (2., 1.9, 0.5, 0.1,
                                               1e-12));
  CHECK (!delay_history_metadata_values_valid (NAN, 1., 0.5, 0.1,
                                                1e-12));
  CHECK (!delay_history_metadata_values_valid (2., INFINITY, 0.5, 0.1,
                                                1e-12));
  CHECK (!delay_history_metadata_values_valid (2., 1., NAN, 0.1,
                                                1e-12));
  CHECK (!delay_history_metadata_values_valid (2., 1., 0.5, INFINITY,
                                                1e-12));
  CHECK (!delay_history_metadata_values_valid (2., 2.1, 0.5, 0.1,
                                                1e-12));
}

int main (void)
{
  test_zero_delay_is_exact();
  test_startup_and_interpolation();
  test_wrapped_ring();
  test_checkpoint_manifests();
  test_variable_timestep_stagger();
  test_metadata_identity_values();
  puts ("delay-history unit tests: PASS");
  return 0;
}
