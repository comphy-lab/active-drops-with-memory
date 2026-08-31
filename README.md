# Active drops with memory

Open Basilisk C implementation of planar chemically fuelled active drops. The
solver couples a CLSVOF interface, incompressible flow, interfacial chemical
production, advection--diffusion and concentration-dependent surface tension.
It supports one, two or four drops, runtime parameter files and checkpoint
restart.

## Requirements

- [Basilisk](https://basilisk.fr/) with `qcc`
- A C compiler and standard maths library
- MPI for parallel builds

The repository includes the local activity model in `src-local/activity.h`.
The benchmark solver uses Basilisk's upstream `two-phase-clsvof.h`;
`src-local/two-phase-clsvof-VP.h` is an experimental viscoplastic extension
and is not included by `dropMove.c`.

## Build

Serial:

```sh
qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
```

MPI:

```sh
CC99=mpicc qcc -O2 -Wall -D_MPI=1 -disable-dimensions \
  dropMove.c -o dropMove -lm
```

## Run

Use a supplied parameter file:

```sh
./dropMove --params parameters/single-drop.cfg
./dropMove --params parameters/two-drop.cfg
./dropMove --params parameters/four-drop.cfg
```

Command-line `key=value` arguments override values loaded earlier:

```sh
./dropMove --params parameters/four-drop.cfg max_level=14 tmax=2 \
  output_dir=output/smoke
```

Supported keys are `Oh`, `Ca`, `Pe`, `AcNum`, `domain`, `max_level`,
`min_level`, `tmax`, `tsnap`, `drops`, `drops_x`, `drops_y`, `drop_radius`,
`drop_spacing`, `ke_limit`, `movement_threshold`, `t0`, `delay_samples`,
`output_dir`, `restart_file` and `resume`. The `drops` shorthand accepts 1, 2
or 4; `drops_x` and `drops_y` define larger rectangular arrays for scale-up
studies. For compatibility with the original Péclet-number scan, a single
positional number is interpreted as `Pe`.

The supplied files demonstrate the one-, two- and four-drop initialisations.
The four-drop file is the large-domain benchmark starting point; its resolution
and domain should be increased together when constructing larger physical
scale-up cases.

## Delayed concentration response

Set `t0` to the response delay. Surface tension then samples the Eulerian
concentration history,

$$
\sigma(\mathbf{x},t) = Ca^{-1} + 4c_L(\mathbf{x},t-t_0).
$$

This is an Eulerian hard-delay model. It is a controlled diagnostic for the
stationary reference problem, not a material-interface memory law: on a moving
or deforming drop, sampling an old field at today's spatial point is not the
same as following an interfacial parcel.

`t0=0` is the instantaneous model: it allocates no history fields and follows
the original surface-tension and adaptation paths exactly. For `t0>0`,
`delay_samples` (default 16, allowed range 2--256) sets the number of nominal
sampling intervals across one delay. The solver keeps a fixed ring of
`delay_samples + 2` adaptive concentration fields and linearly interpolates
between the two stored tracer-time levels surrounding `t-t0`. The recorded
times include Basilisk's half-step tracer staggering and need not be perfectly
uniform when the CFL timestep changes. Before enough history exists, the
prehistory is defined as `c_L(t<0)=c_L(0)`.

All live history fields participate in the concentration AMR criterion. This
prevents a past concentration structure from being coarsened away before it
is sampled, at the cost of additional memory and potentially more refined
cells. Results with positive delay should therefore include a temporal-history
convergence check in `delay_samples` as well as the usual grid study.

## Restart

Each output event writes a rolling checkpoint at
`<output_dir>/restart`. Resume without deleting previous snapshots using:

```sh
./dropMove --params parameters/four-drop.cfg resume=1
```

Set `restart_file=/path/to/checkpoint` to resume from another checkpoint.
Only MPI rank zero creates the output directory; all ranks participate in
parallel `dump()` and `restore()`. A resumed job must use the same physical
parameters as the job that wrote the checkpoint.

For `t0>0`, every dump has a companion `<dump>.delay-history` file containing
the ring ordering, sample times and tracer time represented by the checkpoint.
Both files and the complete named history-field set are required for restart,
and the resumed run must use the same `t0` and `delay_samples`; a mismatch
fails closed. Instantaneous (`t0=0`) checkpoints retain the original
single-file restart behaviour and remove any stale delay sidecar.

## Software test

The deterministic ring ordering, interpolation, checkpoint-manifest,
variable-timestep staggering and metadata contracts can be checked without
running a simulation:

```sh
cc -std=c99 -Wall -Wextra -Werror -pedantic \
  tests/test_delay_history.c -o /tmp/test-delay-history -lm
/tmp/test-delay-history
```

## Output and post-processing

Snapshots and `log.dat` are written below `output_dir`. The
`postProcess-contour/` and `postProcess-vectors/` directories contain source
and Python scripts for extracting interface and velocity data. Generated
binaries, simulation output and macOS metadata are ignored.

## Licence

GPL-3.0. See [LICENSE](LICENSE).
