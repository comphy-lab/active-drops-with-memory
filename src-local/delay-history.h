/**
# Fixed-delay history clock

Small, Basilisk-independent helpers for a bounded ring of uniformly sampled
field histories. The production solver stores one adaptive scalar field per
slot; this header owns only the temporal ordering and interpolation weights.
*/

#ifndef ACTIVE_DROPS_DELAY_HISTORY_H
#define ACTIVE_DROPS_DELAY_HISTORY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double * times;
  int capacity;
  int count;
  int next;
} DelayHistoryClock;

typedef struct {
  int lower;
  int upper;
  double weight;
} DelayHistoryBracket;

static inline void delay_history_clock_reset (DelayHistoryClock * clock)
{
  clock->count = 0;
  clock->next = 0;
  for (int slot = 0; slot < clock->capacity; slot++)
    clock->times[slot] = NAN;
}

static inline int delay_history_oldest (const DelayHistoryClock * clock)
{
  return clock->count == clock->capacity ? clock->next : 0;
}

static inline int delay_history_newest (const DelayHistoryClock * clock)
{
  return clock->count ?
    (clock->next + clock->capacity - 1)%clock->capacity : -1;
}

static inline int delay_history_record (DelayHistoryClock * clock,
                                        double sample_time)
{
  if (!clock || !clock->times || clock->capacity < 1 ||
      !isfinite(sample_time))
    return -1;
  int slot = clock->next;
  clock->times[slot] = sample_time;
  clock->next = (slot + 1)%clock->capacity;
  if (clock->count < clock->capacity)
    clock->count++;
  return slot;
}

/**
Find the two stored times bracketing `target`. Targets before the first
sample use a constant prehistory equal to that first sample. A target newer
than the latest sample is rejected rather than silently shortening the delay.
*/
static inline int delay_history_bracket (const DelayHistoryClock * clock,
                                         double target, double tolerance,
                                         DelayHistoryBracket * bracket)
{
  if (!clock || !clock->times || !bracket || clock->count < 1 ||
      clock->capacity < clock->count || !isfinite(target) ||
      tolerance < 0.)
    return 0;

  int oldest = delay_history_oldest (clock);
  double first = clock->times[oldest];
  if (!isfinite(first))
    return 0;
  if (target < first - tolerance && fabs(first) > tolerance)
    return 0;
  if (target <= first + tolerance) {
    bracket->lower = bracket->upper = oldest;
    bracket->weight = 0.;
    return 1;
  }

  int previous = oldest;
  for (int offset = 1; offset < clock->count; offset++) {
    int current = (oldest + offset)%clock->capacity;
    double lower_time = clock->times[previous];
    double upper_time = clock->times[current];
    if (!isfinite(upper_time) || upper_time <= lower_time)
      return 0;
    if (target <= upper_time + tolerance) {
      if (fabs(target - upper_time) <= tolerance) {
        bracket->lower = bracket->upper = current;
        bracket->weight = 0.;
      }
      else {
        bracket->lower = previous;
        bracket->upper = current;
        bracket->weight = (target - lower_time)/(upper_time - lower_time);
        if (bracket->weight < 0.)
          bracket->weight = 0.;
        else if (bracket->weight > 1.)
          bracket->weight = 1.;
      }
      return 1;
    }
    previous = current;
  }

  int newest = delay_history_newest (clock);
  if (newest >= 0 && fabs(target - clock->times[newest]) <= tolerance) {
    bracket->lower = bracket->upper = newest;
    bracket->weight = 0.;
    return 1;
  }
  return 0;
}

static inline double delay_history_value (double current,
                                          double lower, double upper,
                                          double weight, double delay_time)
{
  return delay_time == 0. ? current : lower + weight*(upper - lower);
}

static inline int delay_history_metadata_values_valid (
  double checkpoint_time, double checkpoint_concentration_time,
  double delay_time, double sample_interval, double tolerance)
{
  return isfinite(checkpoint_time) && checkpoint_time >= 0. &&
    isfinite(checkpoint_concentration_time) &&
    checkpoint_concentration_time >= 0. &&
    isfinite(tolerance) && tolerance >= 0. &&
    checkpoint_concentration_time <= checkpoint_time + tolerance &&
    isfinite(delay_time) && delay_time > 0. &&
    isfinite(sample_interval) && sample_interval > 0.;
}

/**
Validate the named-field portion of a Basilisk dump manifest. The stream must
be positioned immediately after the dump header and `field_count` includes
the first subtree-size field. History fields must be either absent
(`expected_history == 0`) or exactly `cL_history_0000...` with no duplicates.
*/
static inline int delay_history_manifest_matches (FILE * input,
                                                  long field_count,
                                                  int expected_history)
{
  if (!input || field_count < 1 || field_count > 100000 ||
      expected_history < 0)
    return 0;
  int * seen = expected_history ?
    (int *) calloc (expected_history, sizeof(int)) : NULL;
  if (expected_history && !seen)
    return 0;
  int found = 0, ok = 1;
  for (long field_number = 0; ok && field_number < field_count;
       field_number++) {
    unsigned length = 0;
    if (fread (&length, sizeof(length), 1, input) != 1 ||
        length < 1 || length > 1024) {
      ok = 0;
      break;
    }
    char * name = (char *) malloc ((size_t) length + 1);
    if (!name || fread (name, sizeof(char), length, input) != length) {
      free (name);
      ok = 0;
      break;
    }
    name[length] = '\0';
    if (field_number > 0 && !strncmp (name, "cL_history_", 11)) {
      int matched = -1;
      for (int slot = 0; slot < expected_history; slot++) {
        char wanted[64];
        snprintf (wanted, sizeof(wanted), "cL_history_%04d", slot);
        if (!strcmp (name, wanted)) {
          matched = slot;
          break;
        }
      }
      if (matched < 0 || seen[matched])
        ok = 0;
      else {
        seen[matched] = 1;
        found++;
      }
    }
    free (name);
  }
  free (seen);
  return ok && found == expected_history;
}

#endif
