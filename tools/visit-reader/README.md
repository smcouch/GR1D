# GR1D VisIt Restart Reader

This directory contains a VisIt database plugin for GR1D restart HDF5 files.
It is modeled on the BANGR reader under `~/Code/BANGR/tools/visit-reader`,
but it reads GR1D's flat restart layout directly.

## Scope

- Reads `restart_nt_*_time_*.h5` and `restart_*.h5`.
- Treats each restart file as one timestep in an MTMD database.
- Uses `/time` for VisIt time and `/nt` for cycle.
- Uses `/x1` as the radial coordinate.
- Hides ghost zones by default by counting leading negative-radius `x1` cells
  and cropping the same number of outer cells.
- Exposes rank-1 datasets with length `/n1` as VisIt curve variables.
- Exposes M1 arrays as indexed curve variables, for example
  `q_M1_s1_g01_m1`, `q_M1_fluid_s2_g18_m3`, and `eas_s3_g18_k1`.

## Build

On this machine the VisIt XML generators fail with a Qt processor-feature error,
so the generated registration files are checked in manually.

```bash
cd tools/visit-reader
mkdir -p build
cd build
cmake ..
make -j
```

The CMake file installs the plugin to:

```text
~/.visit/3.4.2/darwin-arm64/plugins/databases/
```

It also builds `gr1d_visit_reader_stub`, a small layout checker that does not
require launching VisIt.

```bash
./gr1d_visit_reader_stub ../../../test/restart_nt_0006333247_time_0.9150000.h5
```

## Usage

Open a GR1D restart file in VisIt and select the `GR1DRestartHDF5` database
plugin if it is not auto-detected. Add a Curve plot for variables such as
`rho`, `press`, `temperature`, `v1`, `v_turb`, or an indexed M1 variable.
