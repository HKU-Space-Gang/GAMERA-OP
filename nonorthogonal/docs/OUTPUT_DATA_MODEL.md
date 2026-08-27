# Yin-Yang output data model

The production model does not store one rectangular Cartesian volume. It
stores two overlapping spherical patches plus one coupled ionosphere file.
Consumers must preserve patch identity and use the overlap rule; concatenating
p0 and p1 along an array axis is incorrect.

The supported MHD and automatically derived MI dimensions are summarized in
[Grid and MPI options](GRID_OPTIONS.md).

## Patch topology and coordinates

- `p0` is the Yin patch and `p1` is the Yang patch.
- Each patch has logical dimensions `(NI, NJ, NK)` ordered as radial,
  colatitude-like and longitude-like indices.
- The accepted production dimensions are `(128, 48, 144)` for each patch.
- Both patches cover the full radial interval. Their angular regions overlap
  so that the pair covers the sphere without a polar singularity.
- Grid HDF5 files store physical Cartesian coordinates, not just logical
  angles. With the supplied normalization, one coordinate unit is one Earth
  radius.
- For a physical point `(x,y,z)`, the reference reader uses `(x,y,z)` in the
  Yin frame and `(-x,z,y)` in the Yang frame before forming Yang spherical
  coordinates. It selects the valid patch with the larger angular interior
  margin in the overlap.

The model-native Earth-to-Sun vector is `(-1,0,0)`. The standard figures plot
`SM-X=-x_native`, so the displayed dayside/upstream is at positive SM-X and
the tail is at negative SM-X.

## Compact MHD grid files

There is one static file per patch:

```text
analysis_grid_p0.h5
analysis_grid_p1.h5
```

Their principal datasets are:

| Dataset | Shape | Meaning |
| --- | --- | --- |
| `coordinates` | `(3,NI,NJ,NK)` | Cartesian cell centers |
| `vertices` | `(3,NI+1,NJ+1,NK+1)` | Cartesian cell vertices |
| `global_cells` | `(3,)` | `NI,NJ,NK` |
| `logical_lower`, `logical_upper` | `(3,)` | logical patch bounds |
| `patch_id` | `(1,)` | 0 Yin or 1 Yang |
| `analysis_schema` | `(1,)` | compact MHD schema version |

The first coordinate index is Cartesian component `(x,y,z)`; it is not an
MPI-rank dimension.

## Compact MHD state files

Each physical output time has a matched pair:

```text
analysis_p0_NNNNNN.h5
analysis_p1_NNNNNN.h5
```

Important datasets are:

| Dataset | Shape | Meaning |
| --- | --- | --- |
| `state` | `(8,NI,NJ,NK)` | cell-centered primitive fields |
| `field_names` | `(8,)` | `rho,p,vx,vy,vz,Bx,By,Bz` |
| `time_code` | `(1,)` | normalized solver time |
| `time_seconds` | `(1,)` | physical seconds |
| `normalization` | `(10,)` | `mu0,x,u,time,rho,p,B,moment,omega,sigma` scales |

Never assume a field index without reading `field_names`. Values in `state`
are normalized code units. Multiply density, pressure, velocity and magnetic
field by `normalization[4]`, `[5]`, `[2]` and `[6]`, respectively, to obtain
SI values. The standard renderer then converts to cm^-3, nPa, km/s and nT.

These files are collective, globally assembled patch arrays written with
parallel HDF5. They are different from `restart/checkpoint_*/restart_p*.h5`,
which are rank-local solver restart states and are not plotting inputs.

## Schema-3 MI files

Each `mi_ionosphere_NNNNNN.h5` contains both hemispheres. The shared root
coordinate arrays are:

| Dataset | Shape | Meaning |
| --- | --- | --- |
| `longitude_rad` | `(Nphi,)` | common North/South longitude axis |
| `colatitude_rad` | `(Ntheta,)` | pole-to-low-latitude axis |

The accepted 128 run uses `Nphi=384`, `Ntheta=43`. Each `north/` and `south/`
group contains `(Ntheta,Nphi)` arrays:

| Dataset | Unit | Meaning |
| --- | --- | --- |
| `potential_V` | V | electrostatic potential |
| `fac_parallel_A_m2` | A m^-2 | stored parallel FAC |
| `Sigma_P_S`, `Sigma_H_S` | S | total Pedersen and Hall conductance |
| `Sigma_P_EUV_S`, `Sigma_H_EUV_S` | S | EUV-only contributions |
| `F_E_erg_cm2_s` | erg cm^-2 s^-1 | precipitating electron energy flux |
| `F_N_cm2_s` | cm^-2 s^-1 | precipitating electron number flux |
| `DPB_mask` | dimensionless | diffuse-precipitation mask |
| `DPB_boundary_MLAT_deg` | degrees, shape `(Nphi,)` | DPB boundary |
| `owner_patch` | integer | contributing Yin/Yang patch |

Diagnostic/audit fields such as `fac_yin_A_m2`, `fac_yang_A_m2`, validity,
margin and weight arrays are also present. Root attributes record
`schema_version=3`, physical time, coupling settings, F10.7, precipitation
beta and `sunward_direction_x/y/z`.

The North and South longitude axes are identical and must not be reflected.
Stored longitude maps as follows:

| Stored longitude | MLT | Standard screen position |
| --- | --- | --- |
| 0 | 00 midnight | bottom |
| pi/2 | 06 dawn | right |
| pi | 12 noon | top |
| 3pi/2 | 18 dusk | left |

For upward-current plots, multiply `fac_parallel_A_m2` by the hemisphere
group attribute `upward_fac_multiplier` (normally -1 North and +1 South).
This sign handling is separate from the longitude coordinate and is not a
reason to mirror a hemisphere.

## Reading rule used by the standard scripts

1. Match p0 and p1 using `time_seconds`, not only the filename sequence.
2. Round physical seconds only for cadence matching and endpoint
   deduplication; retain the original time in output metadata.
3. Read `field_names` and `normalization` from the file.
4. Interpolate both patches to the requested Cartesian plane.
5. In the overlap, use the valid patch farther from its angular boundary.
6. Match the MI file by its root `time_seconds` attribute.
7. Verify every required final dataset is finite before accepting products.

The implementation of these rules is frozen in
`scripts/diagnostics/postprocess_compact_base.py`.
