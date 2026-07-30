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
`drop_spacing`, `ke_limit`, `movement_threshold`, `output_dir`, `restart_file`
and `resume`. The `drops` shorthand accepts 1, 2 or 4; `drops_x` and `drops_y`
define larger rectangular arrays for scale-up studies. For compatibility with
the original Péclet-number scan, a single positional number is interpreted as
`Pe`.

The supplied files demonstrate the one-, two- and four-drop initialisations.
The four-drop file is the large-domain benchmark starting point; its resolution
and domain should be increased together when constructing larger physical
scale-up cases.

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

## Output and post-processing

Snapshots and `log.dat` are written below `output_dir`. The
`postProcess-contour/` and `postProcess-vectors/` directories contain source
and Python scripts for extracting interface and velocity data. Generated
binaries, simulation output and macOS metadata are ignored.

## Licence

GPL-3.0. See [LICENSE](LICENSE).
