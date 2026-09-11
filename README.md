# Active drops with finite exterior fuel

Deformable liquid drops consume a finite supply of fuel in the surrounding
fluid. Local depletion raises surface tension and can generate Marangoni
flow. The interface remains a material fluid interface: exhaustion removes
the chemical supply, not the drop itself. Residual flow can relax after
uptake becomes small.

This implementation is a continuum analogue of an object using surrounding
nutrients. It does not represent cell growth, death, division, membranes or
metabolism. It also does not identify fuel with a particular surfactant or
micellar species. The current solver contains no delayed concentration
sampling. The historical constant-production implementation remains in
`src-local/activity.h`, which `dropMove.c` does not include.

## Dimensionless model

Let the reference length, speed and fuel concentration be $R_*$, $U_*$ and
$F_*$. Time is measured in $R_*/U_*$, pressure in $\mu_o U_*/R_*$ and
surface tension in $\mu_o U_*$. Subscripts $o$ and $i$ denote the exterior
and drop fluids. The drop volume fraction is $f$; exterior fraction is
$h=1-f$. The evolved fuel inventory is $q=hF$, where $F$ is exterior fuel
availability normalized by the fixed reference $F_*$. Initially
$F=\mathrm{InitialFuel}$ outside every drop and there is no interior fuel.

The one-fluid momentum equations are

$$
\nabla\cdot\mathbf u=0,\qquad
Re\,\rho_r(f)(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u)
=-\nabla p+\nabla\cdot[\mu_r(f)(\nabla\mathbf u+\nabla\mathbf u^T)]
+\nabla\cdot[\gamma(\mathbf I-\mathbf n\mathbf n)\delta_\Gamma].
$$

Here $\rho_r=1+(\rho_i/\rho_o-1)f$ and
$\mu_r=1+(\mu_i/\mu_o-1)f$. Both temporal and convective inertia are
retained. Small `Re` approximates a low-inertia regime; it is not an
instantaneous Stokes solver.

The fuel equation and uptake law are

$$
\partial_t(hF)+\nabla\cdot(\mathbf u hF)
=\frac1{Pe}\nabla\cdot(h\nabla F)-J(F_s)\delta_\Gamma,
\qquad
J(F_s)=\frac{Da}{Pe}\frac{F_s}{K+F_s},
$$

$$
\gamma=\frac1{Ca}+\mathrm{GammaSlope}(1-F_s).
$$

$F_s$ denotes the exterior interfacial value. Uptake vanishes at zero fuel
and saturates at large fuel availability. With the normal directed from
the drop into the exterior, the intended sharp-interface condition is
$Pe^{-1}\partial_nF=J(F_s)$. There is no source, replenishment or fuel
reservoir imposed at the boundaries: both coordinates are periodic.

| Parameter | Definition or role | Default |
| --- | --- | --- |
| `Re` | $\rho_o U_*R_*/\mu_o$ | 0.01 |
| `Ca` | $\mu_o U_*/\gamma_*$, tension at $F=1$ | 0.1 |
| `Pe` | $U_*R_*/D$ | 1.6 |
| `Da` | $J_{\max}R_*/(DF_*)$ | 1 |
| `FuelHalfSaturation` | $K=K_\mathrm{dim}/F_*>0$ | 0.1 |
| `GammaSlope` | Positive magnitude $-(\partial\gamma_\mathrm{dim}/\partial F_\mathrm{dim})F_*/(\mu_o U_*)$ | 4 |
| `InitialFuel` | Initial exterior availability, between 0 and 1 | 1 |
| `viscosityRatio` | $\mu_i/\mu_o$ | 1 |
| `densityRatio` | $\rho_i/\rho_o$ | 1 |

`GammaSlope` is an independently specified coupling magnitude, not a
parameter named “Ma”. Its default of four is consistent with using the
equal-viscosity planar mobility scale; changing viscosity ratio does not
automatically change this coefficient. `Oh` is derived as
$\sqrt{Ca/Re}$ using the reference tension. Legacy `Oh` and `AcNum`
arguments are rejected rather than silently reinterpreted.

## Discretization and limits

The solver uses Basilisk's CLSVOF interface and integral surface-stress
formulation. Normal capillary and tangential Marangoni stresses act on a
deformable drop. Fuel is a permanent inverse VOF tracer. Its phase-aware
refinement and first-order concentration reconstruction conserve inventory
without assigning fuel to the drop interior.

Diffusion uses exterior cell capacity and nonnegative harmonic exterior
face weights. A face is disconnected when either adjoining exterior
fraction is zero. Inactive algebraic rows are anchored at zero and cannot
exchange fuel. These face weights are an approximation to exterior
connectivity, not a geometrically exact cut-interface diffusion operator.

Uptake uses the reconstructed PLIC segment length in each mixed cell and
the local exterior cell-average concentration. Exactly face-aligned
segments between pure cells are assigned once to the exterior cell.
A backward-Euler Michaelis
update is solved with its positive quadratic root in inventory space. The
same decrement is accumulated as consumed fuel. A post-VOF repair handles
only machine-roundoff negative or orphan inventory. It validates the field
before removing those residues and redistributes their signed inventory
through a global rescaling of positive exterior fuel. No physical source
or unrecorded correction supplies fuel.

With machine epsilon $\epsilon$, the per-cell residue limit is
$64\epsilon\,\mathrm{InitialFuel}$. Absolute transferred inventory is
limited to $256\epsilon M_0$ per repair and $4096\epsilon M_0$ cumulatively,
where $M_0$ is the initial fuel inventory. Transfer is recorded separately
from chemical consumption and persisted in checkpoints. Larger negative
fuel, geometric overshoots, or exhausted correction budgets stop the run.
Zero initial fuel has no correction allowance. Conservation is checked
against the floating-point summation bound; diffusion, reaction and
adaptation retain their strict fuel-state checks.

The surface-tension calculation extends an exterior numerator and support
weight together, then takes their ratio in the nearby drop interior. This preserves a
constant exterior field and avoids treating zero interior inventory as
zero interfacial availability. It remains a low-order extension: neither
the surface trace nor the interfacial uptake is claimed to be a converged
sharp-interface treatment. Geometry, timestep, spatial-resolution and
low-inertia studies are required before interpreting propulsion or arrest.

The geometry is planar two-dimensional, with periodic copies of the drop
configuration. It is not an isolated three-dimensional sphere. A symmetric
single-drop initial condition can remain stationary. Finite-time motion
does not establish an instability threshold; the legacy `PeScan.py` does
not provide a qualified threshold search for this depleting system.

## Build and parameters

Use a project-local installation of
[comphy-lab/basilisk-C](https://github.com/comphy-lab/basilisk-C), pinned to
`v2026-08-30`, with its compiler under `basilisk/src` and installation record
in `basilisk/.comphy-lock`. The generated `.project_config` selects it.

```sh
source .project_config
qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove --params parameters/single-drop.cfg
```

`runCases.sh` builds and runs a selected parameter file. The supplied files
define one, two and four drops. They are example configurations, not
validated production settings.

```sh
./runCases.sh parameters/two-drop.cfg
./dropMove --params parameters/four-drop.cfg Pe=2 Da=0.5
```

All parameters use the same `key=value` parser. Later arguments override
earlier file values. In addition to the model parameters above, the parser
accepts `domain`, `max_level`, `min_level`, `tmax`, `tsnap`, `drops`,
`drops_x`, `drops_y`, `drop_radius`, `drop_spacing`, `ke_limit`,
`movement_threshold`, `output_dir`, `restart_file` and `resume`.
The `drops` shorthand accepts 1, 2 or 4. The default radius is one reference
length. `movement_threshold` is disabled by default; its legacy aggregate
centroid diagnostic is not a per-drop tracker and is unsuitable across a
periodic seam.

## Output, accounting and restart

`log.dat` contains iteration, time, kinetic energy, remaining fuel,
cumulative consumed fuel, inventory defect, drop area and cumulative
absolute roundoff transfer. Roundoff transfer is not physical uptake. The inventory
defect is

$$
\int q\,dV+M_\mathrm{consumed}-M_\mathrm{initial}.
$$

Snapshots contain `fuel` ($q$), `availability` ($F$), `surface_fuel`, `f`,
`d`, `sigmaf`, velocity and checkpoint metadata. The two `getData.c`
extractors output `fuel` directly; multiplying it by $1-f$ again would be
incorrect. Plot labels and scales should identify inventory rather than
the historical produced-species concentration.

```sh
./dropMove --params parameters/single-drop.cfg resume=1
```

Every checkpoint embeds the model version, physical parameters, initial
layout, pressure, face velocities and fuel-accounting totals. Legacy checkpoints and mismatched
physics are rejected. Output location and final time can change on resume.
Restart restores the beginning-of-step state after initialization, then
re-arms its pending solver events. The restart contract preserves physical
time advancement, fuel accounting and numerical continuation; it does not
claim bitwise-identical trajectories across machines or parallel layouts.

## Software checks

```sh
python3 -m unittest discover -s testCases -v
```

The operator probes compile the actual case but replace its time-loop
entry. They check initialization, parameter mapping, positivity, inventory
balance, diffusion, VOF transport, adaptation and restart. Michaelis uptake
is compared with an independent bisection solution of its implicit update.
A separate four-step run of the unmodified main checks uninterrupted
against repeated checkpoint continuation, on uniform and adaptive grids.
A 100-step passive VOF translation at directional CFL 0.49 exercises the
bounded redistribution; larger deliberate violations must be rejected.
They do not validate active-drop motion, uptake kinetics or biological
behaviour against experiments.

## Context and related models

[Michelin (2023), *Self-Propulsion of Chemically Active Droplets*, Annual
Review of Fluid Mechanics 55, 77–101](https://doi.org/10.1146/annurev-fluid-120720-012204)
reviews chemically driven droplet motion, including consumption of a
tension-lowering species and production of a tension-raising species.
[Li and Koch (2022), Journal of Fluid Mechanics 952,
A16](https://doi.org/10.1017/jfm.2022.891) provides a related exterior
surfactant-consumption and planar-mobility setting. The saturating,
finite-inventory law here is a stated modelling choice, not a reproduction
of either paper's complete model.

[ActiveDrops](https://github.com/comphy-lab/ActiveDrops) is the companion
classical active-drop codebase.
[stochastic-active-drops](https://github.com/comphy-lab/stochastic-active-drops)
develops the independent stochastic programme. This repository retains its
historical name and the broader solubilization context while this solver
implements the finite exterior-fuel model described above.

GPL-3.0; see [LICENSE](LICENSE).
