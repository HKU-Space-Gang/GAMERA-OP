# Frozen non-orthogonal Earth production solver

This directory is an independent production subtree of GAMERA-OP. It does not
replace or modify the standard orthogonal solver in the repository root.

The release contains the complete non-orthogonal Yin-Yang global
magnetosphere path needed by a user:

- finite-volume MHD in Cartesian physical components on a non-orthogonal
  Yin-Yang grid;
- the production radial map from 2.5 to 200 Earth radii: a linear inner branch
  joined C1 at 14 RE and logical radius 23/72 to one exponential stretch; the
  72 grid has 0.5 RE inner cells and the 128 grid samples that branch at
  0.28125 RE;
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
3. Choose a supported [grid and MPI preset](docs/GRID_OPTIONS.md).
4. Configure an [idealized or observed-event solar-wind input](docs/SOLAR_WIND.md).
5. Read the [Yin-Yang output data model](docs/OUTPUT_DATA_MODEL.md).
6. Learn the [standard diagnostics and movie workflow](docs/DIAGNOSTICS.md).
7. Only then use the [Sugon 128-resolution example](docs/SUGON.md).

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
│   ├── validate_smoke.py
│   └── diagnostics/              frozen standard renderers
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

The supported 32, 72 and 128 presets, their tested decompositions and their
automatically derived MI grids are listed in
[Grid and MPI options](docs/GRID_OPTIONS.md). The 32 grid is an installation
smoke, the 72 grid is a medium validation option, and only the 128 grid is the
frozen production configuration.

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

## Yin-Yang and MI data structure

MHD output is stored as a matched pair, not as one Cartesian array:

```text
analysis_grid_p0.h5              Yin static grid
analysis_grid_p1.h5              Yang static grid
analysis_p0_NNNNNN.h5            Yin compact state at one physical time
analysis_p1_NNNNNN.h5            Yang compact state at the same time
mi_ionosphere_NNNNNN.h5          North and South MI state in one file
```

For each MHD patch, `coordinates` has shape `(3,NI,NJ,NK)`, `vertices` has
shape `(3,NI+1,NJ+1,NK+1)`, and `state` has shape `(8,NI,NJ,NK)`. The eight
fields are named by the file metadata and are
`rho,p,vx,vy,vz,Bx,By,Bz`. p0 and p1 overlap angularly; they must be
interpolated with the Yin-Yang overlap ownership rule and must never be
concatenated along an array dimension.

Each schema-3 MI file contains shared `longitude_rad` and `colatitude_rad`
axes plus `north/` and `south/` groups. Each hemisphere stores potential,
parallel FAC, total and EUV Pedersen/Hall conductance, FE, FN, DPB fields and
patch-ownership audit fields. The North and South longitude axes are
identical: 00 MLT is bottom, 06 right, 12 top and 18 left. Do not mirror the
South data. Apply each hemisphere's `upward_fac_multiplier` only to convert
stored parallel FAC to upward FAC.

See [Yin-Yang output data model](docs/OUTPUT_DATA_MODEL.md) for dataset shapes,
units, normalization indices, native/model coordinates and reading rules.
Use only the [frozen standard diagnostics](docs/DIAGNOSTICS.md) for production
comparisons.

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
hash receipts. The [solar-wind input guide](docs/SOLAR_WIND.md) gives the exact
HDF5 schema, native coordinate convention, idealized-file example, event-time
alignment, propagation settings and validation checklist.

See [production provenance](docs/PRODUCTION_PROVENANCE.md) for the exact
accepted binary and run evidence.
