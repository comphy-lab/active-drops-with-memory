/**
# Advection/diffusion of temperature tracers
The *stracers* list of temperature tracers must be defined by the calling
code. 
copied and inspired from http://basilisk.fr/src/henry.h & http://basilisk.fr/sandbox/farsoiya/phase_change/phase-change.h
*/

attribute {
  scalar phi1, phi2; // private
  double A;  // activity at the interface
  double D; // Diffusivity of tracer in the diffusive phase
}

extern scalar * stracers;
scalar ActivityFlux[];

#include "diffusion.h"

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

event defaults (i = 0)
{
  fprintf (ferr, "Recompiled!!!\n");
  for (scalar s in stracers) {
#if TREE
    s.refine = refine_bilinear;
    set_prolongation (s, refine_bilinear);
    set_restriction (s, restriction_volume_average);
    s.gradient = p.gradient;
#endif // TREE
  }
}

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 3.2 of [Farsoiya et al.,
2021](#farsoiya2021). */

static scalar * phi_tracers = NULL;

event vof (i++)
{
  phi_tracers = f.tracers;
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c f
    $$
    $$
    \phi_2 = c (1-f)
    $$
    */
		  
    foreach() {
      double a = c[];
      phi1[] = a*f[];
      phi2[] = a*(1. - f[]);
    }
  }
}

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;
  for (scalar c in stracers) {
    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */
    
    scalar phi1 = c.phi1, phi2 = c.phi2;
    foreach() {
      c[] = phi1[] + phi2[];
    }
    delete ({phi1, phi2});
    
    scalar volumic_metric[], dirichlet_source_term[];
    face vector diffusion_coef[];

    scalar diracDelta[];
    // foreach(){
    //   if (interfacial(point, f)){
    //     coord n = interface_normal (point, c), p;
    //     double alpha = plane_alpha (c[], n);
    //     double area = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    //     diracDelta[] = area/sq(Delta); //sqrt(sq(fx) + sq(fy));
    //     ActivityFlux[] = diracDelta[]*c.A;
    //   }
    // }

    foreach(){
      // if (interfacial (point, f))
      double fx = (f[1]-f[-1])/(2*Delta);
      double fy = (f[0,1]-f[0,-1])/(2*Delta);
      diracDelta[] = sqrt(sq(fx) + sq(fy));
      ActivityFlux[] = diracDelta[]*c.A;
    }
  
    foreach() {
      volumic_metric[] = cm[];
      dirichlet_source_term[] = cm[]*ActivityFlux[];
    }
    foreach_face(){
      double ff = (f[] + f[-1])/2.;
      double wt = c.inverse ? (1.-ff) : ff;
      diffusion_coef.x[] = fm.x[]*c.D*wt;
    }
    diffusion (c, dt, D = diffusion_coef, r = dirichlet_source_term, theta = volumic_metric);
    }
}
