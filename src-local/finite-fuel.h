/**
# Exterior fuel with finite inventory

The conserved field `fuel` is $(1-f)F$, per total cell volume. Diffusion
uses only connected exterior cells; reaction consumes the local inventory
through a reconstructed interface segment. This planar implementation
does not model micelles, membrane transport or biological metabolism.
*/
#ifndef FINITE_FUEL_H
#define FINITE_FUEL_H

#include "diffusion.h"
#include <float.h>

scalar fuel[], availability[], surface_fuel[];
double fuel_pe = 1.6, fuel_da = 1., fuel_half_saturation = 0.1;
double fuel_consumed = 0., fuel_initial_inventory = 0.;
double fuel_roundoff_inventory = 0., fuel_roundoff_scale = 1.;
bool fuel_ready = false;

static double fuel_inventory()
{
  double total = 0.;
  foreach(reduction(+:total))
    total += fuel[]*dv();
  return total;
}

static bool fuel_recover()
{
  int invalid = 0;
  foreach(reduction(max:invalid)) {
    double h = 1. - f[];
    if (!isfinite(fuel[]) || fuel[] < 0. ||
        (h <= 0. && fuel[] != 0.)) invalid = 1;
    availability[] = h > 0. ? fuel[]/h : 0.;
    if (!isfinite(availability[])) invalid = 1;
  }
  boundary({availability});
  return !invalid;
}

/**
VOF can leave machine-roundoff inventory in a cell whose exterior fraction
has rounded to zero. Validate the entire field before correcting these
residues. Positive exterior inventory is rescaled to the original total;
the absolute amount changed is recorded independently of chemical uptake.
The per-cell, per-step and cumulative budgets are 64, 256 and 4096 machine
epsilons times the corresponding concentration or initial-inventory scale.
No correction is available when the initial fuel supply is exactly zero.
*/
static bool fuel_repair_roundoff()
{
  double original = 0., valid = 0., residue = 0.;
  long count = 0;
  int invalid = 0;
  double cell_budget = 64.*DBL_EPSILON*fuel_roundoff_scale;
  foreach(reduction(+:original) reduction(+:valid) reduction(+:residue)
          reduction(+:count) reduction(max:invalid)) {
    count++;
    double h = 1. - f[], q = fuel[];
    if (!isfinite(f[]) || f[] < -64.*DBL_EPSILON ||
        f[] > 1. + 64.*DBL_EPSILON || !isfinite(q)) invalid = 1;
    original += q*dv();
    if (q < 0. || (h <= 0. && q != 0.)) {
      residue += fabs(q)*dv();
      if (fabs(q) > cell_budget) invalid = 1;
    }
    else valid += q*dv();
  }
  if (invalid || !isfinite(original) || !isfinite(valid)) return false;
  if (residue == 0.) return true;
  double inventory_scale = fuel_initial_inventory;
  if (!(original > 0.) || !(valid > 0.) || !(inventory_scale > 0.) ||
      residue > 256.*DBL_EPSILON*inventory_scale) return false;
  double ratio = original/valid, transferred = 0., predicted = 0.;
  foreach(reduction(+:transferred) reduction(+:predicted)) {
    double q = fuel[], h = 1. - f[];
    double corrected = q < 0. || h <= 0. ? 0. : q*ratio;
    transferred += fabs(corrected - q)*dv();
    predicted += corrected*dv();
  }
  if (!isfinite(ratio) || ratio <= 0. || !isfinite(transferred) ||
      transferred > 256.*DBL_EPSILON*inventory_scale ||
      !isfinite(predicted) || fabs(predicted - original) >
      8.*max(1., (double)count)*DBL_EPSILON*inventory_scale ||
      fuel_roundoff_inventory + transferred > 4096.*DBL_EPSILON*inventory_scale)
    return false;
  foreach() {
    double h = 1. - f[];
    fuel[] = fuel[] < 0. || h <= 0. ? 0. : fuel[]*ratio;
  }
  fuel_roundoff_inventory += transferred;
  /* Summation itself has an O(N epsilon) floating-point error bound. */
  return fuel_recover() && fabs(fuel_inventory() - original) <=
    8.*max(1., (double)count)*DBL_EPSILON*inventory_scale;
}

/**
Extend exterior values through the nearby interior by convex neighbour
averaging. Only this auxiliary field is extended; inventory remains zero
inside the drop. Constant exterior fuel is reproduced exactly. This is a
bounded, low-order extension, not a one-sided high-order surface trace.
*/
static void fuel_surface_extension()
{
  scalar numerator[], weight[], next[], next_weight[];
  foreach() {
    weight[] = max(0., 1. - f[]);
    numerator[] = fuel[];
  }
  for (int pass = 0; pass < 8; pass++) {
    /* Interpolate numerator and weight with identical operators. Their
       ratio preserves a constant even in coarse/fine ghost cells. */
    boundary({numerator, weight});
    foreach() {
      next[] = numerator[];
      next_weight[] = weight[];
      if (weight[] <= 0.) {
        double value = 0., total_weight = 0., count = 0.;
        foreach_neighbor(1)
          if (weight[] > 0.) {
            value += numerator[];
            total_weight += weight[];
            count += 1.;
          }
        if (total_weight > 0.) {
          next[] = value/count;
          next_weight[] = total_weight/count;
        }
      }
    }
    foreach() {
      numerator[] = next[];
      weight[] = next_weight[];
    }
  }
  foreach()
    surface_fuel[] = weight[] > 0. ? numerator[]/weight[] : 0.;
  boundary({surface_fuel});
}

/**
Backward-Euler diffusion with exterior phase capacity. Harmonic phase
weighting disconnects absent-phase cells. Their dummy diagonal anchors a
zero auxiliary unknown and cannot hold or exchange fuel.
*/
static bool fuel_diffuse (double step)
{
  if (!fuel_recover()) return false;
  double before = fuel_inventory();
  scalar capacity[];
  face vector conductance[];
  foreach() {
    double h = max(0., 1. - f[]);
    capacity[] = cm[]*(h > 0. ? h : 1.);
  }
  foreach_face() {
    double h0 = max(0., 1. - f[]), h1 = max(0., 1. - f[-1]);
    conductance.x[] = h0 > 0. && h1 > 0. ?
      fm.x[]*(2.*h0*h1/(h0 + h1))/fuel_pe : 0.;
  }
  mgstats mg = diffusion(availability, step, D = conductance, theta = capacity);
  int invalid = !isfinite(mg.resa) || mg.resa > 1e-8;
  foreach(reduction(max:invalid)) {
    double h = max(0., 1. - f[]);
    if (h > 0. && (!isfinite(availability[]) || availability[] < 0.))
      invalid = 1;
    fuel[] = h > 0. ? h*availability[] : 0.;
  }
  double after = fuel_inventory();
  return !invalid && fabs(after - before) <= 1e-8*max(1., before);
}

/**
Positive root of the backward-Euler Michaelis uptake equation in inventory
space. `exposure` is timestep times interface area per cell volume times
Da/Pe. This form remains well conditioned for small exterior fractions.
*/
static double fuel_after_uptake (double q, double h, double exposure, double half)
{
  if (q == 0. || h == 0. || exposure == 0.) return q;
  double kh = half*h, b = kh + exposure - q;
  double discriminant = hypot(b, 2.*sqrt(kh)*sqrt(q));
  return b >= 0. ? (2.*kh*q)/(discriminant + b) :
    (discriminant - b)/2.;
}

static bool fuel_consume (double step)
{
  face vector aligned_area[];
  foreach_face()
    aligned_area.x[] = (f[] == 0. && f[-1] == 1.) ||
      (f[] == 1. && f[-1] == 0.);
  restriction((scalar *){aligned_area});
  double consumed = 0.;
  int invalid = 0;
  foreach(reduction(+:consumed) reduction(max:invalid)) {
    double h = 1. - f[], old = fuel[];
    double area_density = 0.;
    if (f[] > 0. && f[] < 1.) {
      coord n = interface_normal(point, f), centre;
      double alpha = plane_alpha(f[], n);
      area_density = plane_area_center(n, alpha, &centre)/Delta;
    }
    else if (f[] == 0.) {
      /* A face-aligned interface has no mixed cell. Assign its segment
         to the exterior cell, once, without creating interior fuel. */
      foreach_dimension()
        area_density += (aligned_area.x[] + aligned_area.x[1])/Delta;
    }
    if (area_density > 0. && fuel_da > 0.) {
      double exposure = step*area_density*fuel_da/fuel_pe;
      double next = fuel_after_uptake(old, h, exposure, fuel_half_saturation);
      if (!isfinite(next) || next < 0. || next > old) invalid = 1;
      fuel[] = next;
      consumed += (old - next)*dv();
    }
  }
  fuel_consumed += consumed;
  return !invalid && fuel_recover();
}

#endif
