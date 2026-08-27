# Frozen non-orthogonal Earth production solver

This directory is an independent production subtree of GAMERA-OP. It does not
replace or modify the standard orthogonal solver in the repository root.

The release contains the complete non-orthogonal Yin-Yang global
magnetosphere path needed by a user:

- finite-volume MHD in Cartesian physical components on a non-orthogonal
  Yin-Yang grid;
- the production radial map from 2.5 to 200 Earth radii: exactly 0.5 RE cells
  from 2.5 to 14 RE, followed by one C1 exponential stretch;
- time-dependent HDF5 solar-wind input and upstream propagation;
- dipole/background-field treatment and the accepted MFE Yin-Yang interface;
- sparse Yin-Yang overset communication;
- electrostatic magnetosphere-ionosphere coupling;
- Fedder/Kaiju electron precipitation with beta 0.8724646;
- the accepted Zhang-offset hybrid diffuse-precipitation boundary;
- Robinson auroral conductance;
- Kaiju/LOMPE F10.7 EUV conductance, using the model-native Earth-to-Sun
  direction (-1, 0, 0) and the q'(SZA) twilight tail through 120 degrees;
- an independent 2 S minimum on total Pedersen conductance;
- Hall-aware cached-tensor BiCGStab/SSOR MI potential solves;
- compact MHD analysis output, schema-3 MI output, and rank-local restart.

The unfinished H+ / inner-magnetosphere, ring-current, live-feedback, replay,
and profiling research branches are deliberately absent. They are not needed
to reproduce the frozen production model and should be developed in separate
feature branches.

## Choose one starting point

If this is your first time using the code:

1. Follow [Local setup](docs/LOCAL_SETUP.md).
2. Build the code and run the two-rank smoke example.
3. Read the output and configuration notes below.
4. Only then use the [Sugon 128-resolution example](docs/SUGON.md).

The local smoke run verifies installation and the complete coupling path. Its
32x12x32 grid and 20-second duration are not suitable for scientific analysis.
To remain compatible with a normal serial-HDF5 laptop installation, it checks
rank-local MHD restart output plus schema-3 MI output. The Sugon example enables
parallel compact analysis output and reproduces the accepted production
configuration.

## Directory map

```text
nonorthogonal/
├── CMakeLists.txt
├── src/                         numerical and production driver source
│   ├── nonorthogonal/           reusable non-orthogonal kernels
│   └── nonorthogonal_backend/   MPI, Yin-Yang, MI and restart backend
├── tests/                       production regression tests only
├── scripts/
│   ├── build_local.sh
│   ├── run_local_smoke.sh
│   └── validate_smoke.py
├── examples/earth_magnetosphere/
│   ├── common/                  frozen solar-wind input
│   ├── local_smoke/             two-rank installation test
│   └── sugon_128/               256-rank production build/run templates
└── docs/
```

## The two commands used locally

From this directory:

```bash
./scripts/build_local.sh
./scripts/run_local_smoke.sh
```

A successful final line is:

```text
SMOKE_VALIDATION_PASS ...
```

The scripts refuse to overwrite an existing build or run result. To keep more
than one smoke result, provide another absolute or relative result directory:

```bash
./scripts/run_local_smoke.sh run/local_smoke_second
```

## MPI process count

Every Yin-Yang run has two patches. The required MPI rank count is therefore:

```text
2 * proc_dims_i * proc_dims_j * proc_dims_k
```

The local config uses 2 ranks. The Sugon config uses
`2 * 4 * 4 * 8 = 256` ranks. A mismatch is a configuration error.

## Time units

`time_stop`, `output_interval`, and `restart_interval` in the YAML file are
code units. With the supplied normalization, one code-time unit is
`x_Norm / u_Norm = 63.71 s`. MI cadence fields ending in `_s` are already in
physical seconds. The production value `time_stop=282.5302150368859` is
18,000 physical seconds (5 hours).

## Output files

Important products are:

- `analysis_grid_p0.h5`, `analysis_grid_p1.h5`: static compact grids;
- `analysis_p0_NNNNNN.h5`, `analysis_p1_NNNNNN.h5`: compact MHD states;
- `mi_ionosphere_NNNNNN.h5`: both hemispheres of schema-3 MI output;
- `checkpoint_NNNNNN/`: rank-local restart state and manifest;
- the run log and input-hash receipts produced by the example scripts.

At an exactly reached floating-point endpoint, the MHD output scheduler may
write two raw p0/p1 files whose physical times differ only at roundoff level.
Treat them as one physical time and retain the raw files as provenance. MI
output is not duplicated. Standard postprocessing must select the file nearest
each requested physical cadence rather than count raw filenames blindly.

## Frozen versus configurable

For a comparable production study, keep these compile-time choices fixed:

- `COORD_TYPE=3`, `PROBLEM=22`;
- strict Kaiju magnetic and gas walls ON;
- MFE interface ON;
- sparse Yin-Yang overset ON and profiling OFF;
- H(div) reconciliation OFF;
- radial map version 4, inner radius 2.5 RE, outer radius 200 RE.

Solar-wind history, run duration, output cadence and physically motivated MI
parameters are YAML inputs. Change one physical assumption at a time, record
the changed YAML in the run directory, and keep the automatically generated
hash receipts.

See [production provenance](docs/PRODUCTION_PROVENANCE.md) for the exact
accepted binary and run evidence.
