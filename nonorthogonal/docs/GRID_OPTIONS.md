# Grid and MPI options

Every size in this document is written as `NI x NJ x NK x 2`. `NI`, `NJ`
and `NK` are the radial and two angular cell counts **for one patch**; the
final `x 2` denotes the Yin and Yang patches. It is not an additional array
dimension inside one patch file.

## Supported presets

| Preset | MHD grid per patch | Total MHD cells | MPI decomposition per patch | Total MPI ranks | Automatically derived MI grid per hemisphere | Intended use |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Installation smoke | `32 x 12 x 32` | 24,576 | `1 x 1 x 1` | 2 | `96 x 12` | laptop installation and coupling-path check only |
| Medium validation | `72 x 24 x 72` | 248,832 | `2 x 4 x 4` | 64 | `192 x 22` | short, lower-cost A/B and development validation |
| Frozen production | `128 x 48 x 144` | 1,769,472 | `4 x 4 x 8` | 256 | `384 x 43` | production science and published comparisons |

The MHD-cell total already includes both patches. The MI table uses
`Nphi x Ntheta`; its dimensions are derived from the MHD angular spacing,
the mapped polar-cap extent and `mi_angular_oversampling: 2.0`.

Only the 128 preset is the frozen production configuration. A result obtained
with the smoke or 72 preset must be labelled with its actual grid and must not
be presented as a 128-production result.

## Radial map shared by the presets

All three presets sample the same continuous radial-map version 4:

- inner MHD boundary at 2.5 RE;
- a linear inner branch joined at logical radius `23/72` to 14 RE;
- one C1-continuous exponential branch from 14 to the 200 RE outer boundary.

The mapping is calibrated so that the 72 radial grid has exactly 23 uniform
0.5 RE cells from 2.5 through 14 RE. Changing `NI` samples the same continuous
map more or less densely; it does not preserve a fixed count of 0.5 RE cells:

| `NI` | Width of a cell wholly inside the linear branch | Full linear cells before the transition cell |
| ---: | ---: | ---: |
| 32 | 1.125 RE | 10 |
| 72 | 0.5 RE | 23; the knot lies exactly on a cell boundary |
| 128 | 0.28125 RE | 40 |

For `NI=32` and `NI=128`, one cell straddles the C1 join because `NI*23/72`
is not an integer. The 32 smoke grid is deliberately coarse and is not a
magnetosphere science grid. The 128 grid resolves the inner magnetosphere
more finely than the 0.5 RE design target.

Do not change the compile-time inner radius, outer radius, radial-map version
or stretch constant when making a resolution-only comparison.

## Installation smoke: the local default

The command

```bash
./scripts/run_local_smoke.sh
```

uses `examples/earth_magnetosphere/local_smoke/config.yaml` and launches two
MPI ranks. It runs for 20 physical seconds, starts in the Bz=-5 nT interval,
and validates rank-local restart plus schema-3 MI output. It is intentionally
small so a new student can confirm the installation quickly.

Do not turn this smoke configuration into a long science run merely by
increasing `time_stop`: its outer-magnetosphere resolution is insufficient.

## Medium 72 validation preset

The validated 72-grid decomposition is:

```yaml
ni_global: 72
nj_global: 24
nk_global: 72
proc_dims_i: 2
proc_dims_j: 4
proc_dims_k: 4
mi_angular_oversampling: 2.0
```

It requires:

```text
2 patches * 2 * 4 * 4 = 64 MPI ranks
```

Use it on a workstation or cluster allocation that can launch 64 ranks. Start
from a copy of the Sugon 128 example in a new directory, change the six grid
and decomposition values above, and change the run allocation to 64 tasks.
Keep all physical precipitation, conductance, solar-wind and MI parameters
unchanged for a resolution A/B test. The run receipt must say:

```text
grid=72x24x72x2
mpi_ranks=64
proc_dims_per_patch=2x4x4
mi_grid=192x22_per_hemisphere
```

The frozen full two-panel renderer currently accepts the 128 radial grid
only. Do not silently alter that standard script for a 72 run; a 72 renderer
must be reviewed and identified separately from the frozen 128 products.

## Frozen 128 production preset

The checked-in files under
`examples/earth_magnetosphere/sugon_128/` already contain:

```yaml
ni_global: 128
nj_global: 48
nk_global: 144
proc_dims_i: 4
proc_dims_j: 4
proc_dims_k: 8
mi_angular_oversampling: 2.0
```

The required process count is:

```text
2 patches * 4 * 4 * 8 = 256 MPI ranks
```

This preset samples the v4 inner linear branch at 0.28125 RE per full cell;
it does not contain only 23 inner cells. This increased radial sampling is one
reason it is the production-resolution option.

Use the build, run and diagnostics dependency chain described in
[Sugon production setup](SUGON.md). The frozen standard diagnostics and the
accepted performance evidence correspond to this preset.

## Rules for defining another grid

The code checks these conditions at startup:

1. `NI`, `NJ` and `NK` must each be positive.
2. Each global cell count must be exactly divisible by its matching
   `proc_dims` value.
3. The launched MPI size must equal
   `2 * proc_dims_i * proc_dims_j * proc_dims_k`.
4. Both patches use the same global dimensions and decomposition.
5. MI dimensions should normally remain automatic; omit explicit MI counts
   and keep `mi_angular_oversampling: 2.0` for production-like comparisons.

A new grid is a new numerical configuration. Before using it for science,
run the regression tests, a short finite-output gate, a restart gate and an
A/B comparison against the frozen 128 preset. Record the exact dimensions,
rank decomposition, MI grid and source commit in the run receipt.
