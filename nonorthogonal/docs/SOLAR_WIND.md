# Solar-wind input: idealized schedules and observed events

This guide describes the HDF5 solar-wind input used by the frozen
non-orthogonal Earth model. Read it before replacing the supplied upstream
file. The same reader and propagation path are used for an idealized IMF
schedule and for an observed event; only the construction and quality control
of the HDF5 file differ.

The supplied production example is:

```text
nonorthogonal/examples/earth_magnetosphere/common/
    upstream_zero_then_minus5.h5
```

It contains zero IMF through one hour and then `Bz=-5 nT`. Its corresponding
YAML settings are in
[`sugon_128/config.yaml`](../examples/earth_magnetosphere/sugon_128/config.yaml).

## 1. Model-native coordinates

Do not copy conventional GSM/GSE X components into this model without
checking their signs.

- Earth-to-Sun is model-native `(-1, 0, 0)`.
- The upstream boundary is at negative native X.
- Solar wind enters from native `-X` and flows toward native `+X`.
- Therefore a normal Earth-directed input has **positive native `Vx`**.
- The standard plots display `SM-X = -native-X`, so their dayside appears on
  the positive horizontal side even though the upstream computational boundary
  is at negative native X.

For a vector already expressed in the GSM/GSE axes used by the established
event-preparation workflow, the release input convention is:

```text
native X = -source X
native Y =  source Y
native Z =  source Z
```

Accordingly, that workflow changes the signs of `Vx` and `Bx` and retains
`Vy`, `Vz`, `By`, and `Bz`. Southward IMF remains negative `Bz`. First put all
vector components into the intended scientific frame and record that frame;
do not silently combine components from incompatible frames.

## 2. Required HDF5 schema

Every primary dataset must be a non-empty one-dimensional numeric array with
the same length. The loader reads it as double precision.

| Dataset | Physical-input meaning and unit | Required |
|---|---|---|
| `/T` | sample time in seconds | yes |
| `/D` | proton number density in cm^-3 | yes |
| `/Vx`, `/Vy`, `/Vz` | native velocity, m/s or km/s as selected in YAML | yes |
| `/Bx`, `/By`, `/Bz` | native magnetic field in nT | yes |
| `/P` | **thermal** pressure in nPa | `/P` or `/Temp` |
| `/Temp` | proton temperature in K | allowed only when physical units are used and `/P` is absent |

If both `/P` and `/Temp` are present, `/P` takes precedence. With `/Temp`, the
loader computes proton thermal pressure as `D k_B Temp`. It does not add an
electron or helium pressure model.

The input must also satisfy all of the following:

- `/T` is finite and strictly increasing; duplicate times are invalid;
- all state values are finite;
- `/D` and the resulting thermal pressure are strictly positive;
- `Vx` is nonzero, because the ballistic propagation map is singular at zero
  normal speed.

Optional scalar or length-one datasets are:

| Dataset | Meaning |
|---|---|
| `/ByC` | dimensionless IMF-front coefficient associated with `By` |
| `/BzC` | dimensionless IMF-front coefficient associated with `Bz` |
| `/Bx0` | magnetic offset in nT for a physical-units file |

`/ByC` and `/BzC` must appear together. When their relation is enabled, the
reader replaces every input `Bx` by

```text
Bx = Bx0 + ByC * By + BzC * Bz.
```

Extra groups, datasets, and attributes are allowed and are ignored by the
solver. Use them to retain source URLs, UTC epochs, quality masks, processing
versions, and hashes.

### Code-unit files

`wind_input_units: code` is supported for controlled numerical tests. In that
mode `/T`, `/D`, `/P`, velocity, and magnetic field must all already be
normalized, and `/P` is mandatory. `/Temp` cannot replace `/P`. New Earth
science inputs should normally use `wind_input_units: physical` so their units
remain auditable.

## 3. YAML keys

A physical-units input should declare every relevant option explicitly:

```yaml
wind_file: upstream_input.h5
wind_input_units: physical
wind_velocity_units: km/s
wind_interpolation: linear
wind_reference_x: -200.0
wind_reference_y: 0.0
wind_reference_z: 0.0
wind_time_offset: 0.0
```

The meanings are:

- `wind_file`: absolute path, or a path relative to the directory containing
  the YAML file;
- `wind_input_units`: `physical` or `code`;
- `wind_velocity_units`: `m/s` or `km/s` for a physical file;
- `wind_interpolation`: `linear` or `step`;
- `wind_reference_x/y/z`: location represented by the time series, in model
  code coordinates. With the supplied Earth normalization these numbers are
  in RE;
- `wind_time_offset`: offset added to simulation time before sampling. This is
  in **code-time units, not seconds**.

The supplied normalization gives

```text
Time_Norm = x_Norm / u_Norm = 6371000 / 100000 = 63.71 s.
```

Therefore convert a desired physical offset with

```text
wind_time_offset = desired_offset_seconds / 63.71.
```

Do not enter `3600` to mean one hour; one physical hour is approximately
`56.506043` code-time units. The local smoke configuration deliberately uses
a very large positive offset only to clamp the short installation test to the
final southward-IMF state. It is not an event-timing example.

Optional YAML overrides are:

```yaml
wind_by_coefficient: 0.0
wind_bz_coefficient: 0.0
wind_bx_offset: 0.0
wind_enforce_bx_relation: 0
```

The two coefficients must be supplied together. `wind_enforce_bx_relation`
must be `0` or `1`. If `/ByC` and `/BzC` exist in the file, the relation is
enabled by default; otherwise measured/input `Bx` is retained by default. A
YAML coefficient override enables the relation unless
`wind_enforce_bx_relation: 0` is also set.

For a physical event, prefer putting `ByC`, `BzC`, and `Bx0` in the HDF5 file.
The file's `Bx0` is then interpreted in nT and normalized by the loader. A
YAML `wind_bx_offset` is a code-unit override and is not automatically treated
as nT.

## 4. Idealized solar wind

Use `step` interpolation for an intentionally abrupt state schedule and
`linear` interpolation for a prescribed ramp. Always provide padding before
the first requested state and after the end of the run because sampling
outside the file range is clamped to the first or last record.

The following complete Python program creates a physical-units idealized file
with `Bz=0` for 0--1 h and `Bz=-5 nT` for 1--5 h. Save it as
`make_idealized_wind.py`, activate the Python environment described in
[Local setup](LOCAL_SETUP.md), and run `python3 make_idealized_wind.py`.

```python
from pathlib import Path

import h5py
import numpy as np

output = Path("upstream_zero_then_minus5.h5")

# Seconds relative to simulation start. The first and last records are padding.
T = np.array([-86400.0, 0.0, 3600.0, 86400.0])
n = T.size

fields = {
    "T": T,
    "D": np.full(n, 5.0),                 # cm^-3
    "P": np.full(n, 0.00802858523371),    # thermal pressure, nPa
    "Vx": np.full(n, 400.0),              # native km/s; inflow is +X
    "Vy": np.zeros(n),
    "Vz": np.zeros(n),
    "Bx": np.zeros(n),
    "By": np.zeros(n),
    "Bz": np.array([0.0, 0.0, -5.0, -5.0]),
}

with h5py.File(output, "x") as h5:
    for name, values in fields.items():
        h5.create_dataset(name, data=np.asarray(values, dtype=np.float64))
    h5.attrs["description"] = "0--1 h zero IMF; 1--5 h Bz=-5 nT"
    h5.attrs["coordinate_convention"] = "model native: Sun at -X, inflow +Vx"
    h5.attrs["velocity_unit"] = "km/s"
    h5.attrs["magnetic_unit"] = "nT"
    h5.attrs["density_unit"] = "proton cm^-3"
    h5.attrs["pressure_unit"] = "thermal nPa"

print(output.resolve())
```

Use it from YAML as follows:

```yaml
wind_file: upstream_zero_then_minus5.h5
wind_input_units: physical
wind_velocity_units: km/s
wind_interpolation: step
wind_reference_x: -200.0
wind_reference_y: 0.0
wind_reference_z: 0.0
wind_time_offset: 0.0
```

At an exact sample time, `step` interpolation selects the new sample, so the
above file changes to `Bz=-5 nT` at 3600 s. To model a finite transition,
insert records spanning the ramp and use `linear` interpolation.

To change the idealized solar wind, edit the arrays rather than the solver.
For a controlled A/B experiment, change only one physical input at a time and
keep the original HDF5 file and its SHA-256 hash with the result.

## 5. Observed-event solar wind

An event file needs more care than an idealized schedule. Recommended sources
include one-minute [NASA High Resolution OMNI](https://omniweb.gsfc.nasa.gov/ow_min.html)
data for a terrestrial event or calibrated ACE/DSCOVR measurements when a
separate propagation analysis is intended.

### 5.1 Decide what time and position the source represents

Before processing any values, record:

1. data product and version;
2. native cadence and time standard;
3. coordinate frame for every vector;
4. whether timestamps are raw spacecraft times or already propagated target
   arrival times;
5. position represented by the timestamps.

High-resolution OMNI timestamps are already time shifted to an estimated
arrival at the terrestrial bow-shock nose; NASA's
[HRO construction and time-shifting documentation](https://omniweb.gsfc.nasa.gov/html/sc_merge_data1.html)
states that its record times are target-arrival rather than observation times.
Do not treat them as raw L1 times and propagate them a second time. For an
already propagated product, choose a reference position consistent with its
target convention, commonly near native `(0,0,0)`, and document the choice.
For raw spacecraft data, use a physically consistent spacecraft/reference
position and perform or retain the intended propagation exactly once.

### 5.2 Select an event epoch

Store `/T` as seconds relative to a named UTC event epoch:

```text
T = sample_UTC - event_zero_UTC.
```

If solver time zero is not the event epoch, use

```text
wind_time_offset =
    (simulation_start_UTC - event_zero_UTC).total_seconds() / 63.71.
```

Write both UTC values and the formula result into HDF5 attributes and the run
receipt. Remember that a nonzero reference position and tilted front add a
position-dependent ballistic delay, so the file must extend beyond the run on
both sides.

### 5.3 Map and clean the variables

For each record:

1. replace documented source fill values with missing values;
2. convert vectors into one documented scientific frame;
3. apply the model-native X convention described in section 1;
4. map proton density to `/D`;
5. map proton temperature to `/Temp`, or compute a scientifically justified
   **thermal** pressure for `/P`;
6. remove duplicate timestamps and sort by increasing time;
7. handle short gaps with a declared method and retain a quality mask;
8. reject or explicitly segment long gaps rather than silently bridging them;
9. add time padding required by the run and propagation delay.

OMNI's commonly named `Flow Pressure` variable is solar-wind **dynamic/ram
pressure**, as shown by NASA's
[derived-parameter definition](https://omniweb.gsfc.nasa.gov/ftpbrowser/bow_derivation.html).
It must not be written to `/P`. The MHD reader interprets `/P` as thermal
pressure. For OMNI, the safer default is to write proton `/D` and `/Temp` and
omit `/P`, allowing the loader to compute `n k_B T`.

If the event workflow fits the divergence-compatible Lyon relation, store the
fit as `/ByC`, `/BzC`, and `/Bx0` and retain the measured `/Bx` plus fit report
for provenance. Otherwise omit the coefficients or set
`wind_enforce_bx_relation: 0` to preserve measured `Bx`.

### 5.4 Minimal event HDF5 writer

After source-specific fill-value removal, frame conversion, and gap handling,
export a clean CSV with this header:

```text
utc,D_cm3,Temp_K,Vx_source_kms,Vy_source_kms,Vz_source_kms,Bx_source_nT,By_source_nT,Bz_source_nT
```

The following writer converts the conventional source-X sign to the release's
native convention and refuses non-finite, non-positive, duplicate, reversed,
or longer-than-five-minute-gap input. Change `event_zero_utc`, the declared
source frame, and the gap threshold for the actual study. It intentionally
does not fetch data or conceal missing intervals; source-specific preprocessing
and its quality masks remain part of the event provenance.

```python
import csv
from datetime import datetime, timezone

import h5py
import numpy as np

input_csv = "processed_event.csv"
output_h5 = "event_wind.h5"
event_zero_utc = datetime.fromisoformat("2015-03-17T00:00:00+00:00")
max_gap_seconds = 300.0

columns = {
    "D": [], "Temp": [],
    "Vx": [], "Vy": [], "Vz": [],
    "Bx": [], "By": [], "Bz": [],
}
times = []

with open(input_csv, newline="", encoding="utf-8") as stream:
    for row in csv.DictReader(stream):
        utc = datetime.fromisoformat(row["utc"].replace("Z", "+00:00"))
        if utc.tzinfo is None:
            raise ValueError("Every UTC value must include a timezone")
        utc = utc.astimezone(timezone.utc)
        times.append((utc - event_zero_utc).total_seconds())
        columns["D"].append(float(row["D_cm3"]))
        columns["Temp"].append(float(row["Temp_K"]))
        # Established release mapping: native X=-source X; Y/Z unchanged.
        columns["Vx"].append(-float(row["Vx_source_kms"]))
        columns["Vy"].append(float(row["Vy_source_kms"]))
        columns["Vz"].append(float(row["Vz_source_kms"]))
        columns["Bx"].append(-float(row["Bx_source_nT"]))
        columns["By"].append(float(row["By_source_nT"]))
        columns["Bz"].append(float(row["Bz_source_nT"]))

T = np.asarray(times, dtype=np.float64)
data = {name: np.asarray(values, dtype=np.float64)
        for name, values in columns.items()}

if T.size == 0 or any(values.shape != T.shape for values in data.values()):
    raise ValueError("Empty or unequal-length event arrays")
if not np.all(np.isfinite(T)) or not all(np.all(np.isfinite(values))
                                                for values in data.values()):
    raise ValueError("Non-finite event value")
if not np.all(np.diff(T) > 0.0):
    raise ValueError("UTC samples must be strictly increasing and unique")
if T.size > 1 and np.max(np.diff(T)) > max_gap_seconds:
    raise ValueError("Event contains a gap longer than the accepted threshold")
if not np.all(data["D"] > 0.0) or not np.all(data["Temp"] > 0.0):
    raise ValueError("Density and temperature must be positive")
if not np.all(data["Vx"] > 0.0):
    raise ValueError("Expected positive model-native solar-wind Vx")

with h5py.File(output_h5, "x") as h5:
    h5.create_dataset("T", data=T)
    for name, values in data.items():
        h5.create_dataset(name, data=values)
    h5.attrs["event_zero_utc"] = event_zero_utc.isoformat()
    h5.attrs["source_frame"] = "REPLACE_WITH_VERIFIED_SOURCE_FRAME"
    h5.attrs["coordinate_mapping"] = "native X=-source X; Y/Z unchanged"
    h5.attrs["input_csv"] = input_csv
    h5.attrs["velocity_unit"] = "km/s"
    h5.attrs["magnetic_unit"] = "nT"
    h5.attrs["density_unit"] = "proton cm^-3"
    h5.attrs["temperature_unit"] = "K"

print(output_h5)
```

For a real production input, also store or archive the raw-data checksum,
retrieval URL/time, product version, fill-value rules, interpolation method,
quality masks, and the checksum of this generated HDF5 file. If the source CSV
is already in model-native coordinates, remove the two X sign changes and
record that fact instead of applying the mapping twice.

### 5.5 Event YAML example

The following assumes that `/T` is relative to the chosen event epoch, the
data are in physical units, velocities are km/s, and the source has been
prepared at the stated model-native reference position:

```yaml
wind_file: st_patricks_20150317_omni_1min.h5
wind_input_units: physical
wind_velocity_units: km/s
wind_interpolation: linear
wind_reference_x: 0.0
wind_reference_y: 0.0
wind_reference_z: 0.0
wind_time_offset: -70.63255375922147
wind_enforce_bx_relation: 0
```

Here the example offset is `-4500 s / 63.71 s`; it is illustrative and must
be recomputed from the actual event and simulation epochs. Use `linear`
interpolation for normal time-series observations. Never copy the numeric
offset or reference point to another event without re-deriving them.

## 6. Propagation performed by the model

The reader first samples monitor time

```text
t_monitor = t_simulation + wind_time_offset.
```

For `ay=ByC`, `az=BzC`, and the current native `Vx`, the front velocity is

```text
v_front = Vx / (1 + ay^2 + az^2) * (1, -ay, -az).
```

At model point `x` relative to `wind_reference`, the sample time is

```text
delay    = dot(x - reference, v_front) / dot(v_front, v_front)
t_sample = t_monitor - delay.
```

All quantities in these equations are in normalized code units after file
loading. The model clamps `t_sample` outside the input range to the nearest
endpoint. Clamping avoids a crash but is not acceptable as accidental event
coverage; use sufficient input padding and audit the actual delay range.

## 7. Pre-run validation

Inspect the structure and values before submitting an expensive run:

```bash
h5dump -n upstream_input.h5
h5dump -d /T -d /D -d /Vx -d /Bx -d /By -d /Bz upstream_input.h5
sha256sum upstream_input.h5
```

On macOS, `shasum -a 256 upstream_input.h5` is equivalent. Also make a QA
plot of every primary variable versus UTC and verify:

- one-dimensional equal-length arrays;
- finite values and strictly increasing time;
- positive density and thermal pressure or temperature;
- expected native `Vx` sign and no zero `Vx`;
- expected IMF signs after coordinate conversion;
- no source fill values or unexplained discontinuities;
- gap masks and interpolation policy are visible;
- coverage includes the complete run plus propagation padding;
- YAML path, units, interpolation, reference, and code-time offset are correct;
- coefficient/relation choice is explicit;
- raw-source and generated-file hashes are recorded.

Then run a short smoke/gate with the same HDF5 and YAML physics settings. The
log must report the loaded sample range, interpolation mode, front normal, and
whether the Bx relation is enabled. Inspect the first and last boundary states
before authorizing a long production run.

## 8. Common mistakes

- **Negative native `Vx`:** copied from conventional GSM/GSE without the
  release's X mapping.
- **Dynamic pressure in `/P`:** produces the wrong thermal state; use thermal
  pressure or `/Temp`.
- **Offset entered in seconds:** `wind_time_offset` is code time; divide
  physical seconds by `63.71` for the supplied normalization.
- **Double propagation:** using already propagated OMNI timestamps as though
  they were raw L1 observations.
- **Hidden endpoint clamping:** a run appears to work but spends part of its
  interval at a frozen first or last sample.
- **Wrong interpolation:** `step` is appropriate for a designed discontinuous
  schedule, while observed event data normally require `linear`.
- **Partial coefficient override:** `wind_by_coefficient` and
  `wind_bz_coefficient` must be provided together.
- **Untracked preprocessing:** an event HDF5 without source, frame, epoch, gap,
  and hash metadata is not reproducible.
