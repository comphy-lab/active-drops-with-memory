# Finite exterior-fuel active drops

Read `README.md` for the model and dimensionless definitions. This branch
describes a deformable drop consuming surrounding fuel; the delayed-response
implementation is a separate workstream.

## Structure and execution

- `dropMove.c` owns the case, one parameter parser, initialization, diagnostics
  and restart contract. Keep simulation logic visible in the standalone case.
- `src-local/finite-fuel.h` owns exterior inventory transport, uptake and the
  auxiliary surface-field extension. The old `activity.h` is not included.
- `parameters/` contains small one-, two- and four-drop examples.
- `testCases/test_finite_fuel.py` exercises operators and bounded actual
  timesteps; it is not a propulsion or biological validation campaign.

Use the project-local Basilisk v2026-08-30 selected by `.project_config`.
Keep the compiler, generated binaries and simulation output ignored.
Run `python3 -m unittest discover -s testCases -v` after numerical changes.

## Numerical contracts

Fuel inventory is `q=(1-f)F`, with zero interior inventory. Physical uptake
and conservative machine-roundoff redistribution have separate accounting.
Never repair a substantial conservation, positivity or geometry failure by
clipping fields or relaxing a tolerance without evidence.

`Re`, `Ca`, `Pe`, `Da`, `GammaSlope` and material ratios have explicit scale
definitions. `Oh` is derived. Do not reinstate unexplained coupling constants
or confuse reference-velocity Pe with a measured swimming-speed threshold.

The initialization order and saved face velocity, pressure, fuel and scheduler
state are part of the restart contract. Test real uninterrupted-versus-resumed
evolution after changing that contract; field readback alone is insufficient.

This is a planar finite-resource model with approximate interface transport.
It does not model cell growth, division, biological death or delayed sensing.
Preserve the independent solubilizing General and reduced stochastic models.
