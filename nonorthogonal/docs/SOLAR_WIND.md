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

The supported event path is deliberately a single command. Do not download
individual columns through the OMNIWeb form, merge tables, calculate an L1
delay, or hand-edit timestamps. The checked-in tool requests one definitive
NASA file, converts it, validates it, and writes every run-side setting.

### 5.1 NASA product and variables

The tool uses NASA CDAWeb/SPDF HAPI dataset
[`OMNI_HRO2_1MIN`](https://cdaweb.gsfc.nasa.gov/hapi/info?id=OMNI_HRO2_1MIN),
the definitive one-minute high-resolution OMNI product (DOI
[`10.48322/mj0k-fq60`](https://doi.org/10.48322/mj0k-fq60)). Its timestamps
have already been shifted to estimated arrival at the terrestrial bow-shock
nose, as described in NASA's
[high-resolution OMNI documentation](https://omniweb.gsfc.nasa.gov/html/HROdocum.html).
Therefore **do not apply another ACE/DSCOVR L1-to-Earth shift**.

The script requests exactly these NASA variables:

| NASA HAPI field | Frame/unit | GAMERA-OP output | Automatic conversion |
|---|---|---|---|
| `BX_GSE` | GSE, nT | `/Bx` | reverse X sign |
| `BY_GSM`, `BZ_GSM` | GSM, nT | `/By`, `/Bz` | retain Y/Z signs |
| `Vx`, `Vy`, `Vz` | GSE, km/s | `/Vx`, `/Vy`, `/Vz` | reverse X sign only |
| `proton_density` | proton cm^-3 | `/D` | direct |
| `T` | K | `/Temp` | direct |

GSE and GSM share the X axis, which makes the mixed magnetic selection above
well defined. NASA fill values from the HAPI metadata are changed to missing
values before interpolation. The script does not request OMNI `Pressure`:
that field is dynamic/ram pressure, whereas solver `/P` means **thermal**
pressure. Writing `/D` and `/Temp` lets the loader compute proton thermal
pressure consistently.

### 5.2 Install the small Python toolchain once

From `GAMERA-OP/nonorthogonal`, activate the environment created by
[Local setup](LOCAL_SETUP.md), or make one now:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy h5py matplotlib
```

The downloader uses only Python's standard HTTPS client; no NASA account,
browser export, `cdasws`, or SpacePy installation is required.

### 5.3 Complete event example: St. Patrick's Day storm interval

Run this one command from `GAMERA-OP/nonorthogonal`:

```bash
python3 scripts/solar_wind/fetch_nasa_omni_event.py \
  --simulation-start 2015-03-17T20:10:00Z \
  --duration-hours 3 \
  --padding-minutes 60 \
  --output run_inputs/st_patricks_2015/st_patricks_20150317_omni_1min.h5
```

That invocation performs the full workflow:

1. requests `OMNI_HRO2_1MIN` once from NASA CDAWeb HAPI;
2. adds the requested 60 minutes before and after the focused three-hour run;
   for other events the default is 120 minutes;
3. removes the documented source fill values;
4. linearly fills missing runs up to the default 10-minute limit and rejects
   it if any required field has a longer missing run;
5. applies the release's native X-axis convention;
6. writes `/T` in seconds relative to the requested simulation start;
7. validates cadence, finiteness, density, temperature, native inflow, and
   that padding exceeds the worst observed advection delay over 200 RE;
8. creates the HDF5, run settings, QA plot, and checksummed provenance receipt.

On success the last block begins with `NASA_OMNI_EVENT_READY`. The directory
contains:

```text
run_inputs/st_patricks_2015/
├── st_patricks_20150317_omni_1min.h5          solver input
├── st_patricks_20150317_omni_1min.wind.yaml   copy-ready YAML keys
├── st_patricks_20150317_omni_1min.qa.png      converted-variable QA
├── st_patricks_20150317_omni_1min.receipt.json
└── raw/
    ├── NASA_OMNI_HRO2_1MIN_..._UTC.csv        exact NASA response
    └── NASA_OMNI_HRO2_1MIN_..._UTC.csv.source.json
```

This focused interval is used because the public one-minute product contains
a long data outage earlier on 17 March; the tool rejects a full-day input
instead of silently interpolating across it. The selected command has been
run against NASA data: the largest required-field gap is 5 minutes, and the
60-minute padding passes the observed-speed audit for the 200-RE domain.

The raw CSV is the single download. Re-running the same time range uses that
cache; add `--refresh` only when a deliberate new NASA retrieval is wanted.
The receipt records the exact query URL, retrieval time, raw and generated
SHA-256 hashes, maximum gap by variable, and interpolation count.

### 5.4 Run configuration: use the generated file, do not derive a shift

Open the generated `.wind.yaml` and copy its keys into the science
`config.yaml`. For this example it contains:

```yaml
wind_file: st_patricks_20150317_omni_1min.h5
wind_input_units: physical
wind_velocity_units: km/s
wind_interpolation: linear
wind_reference_x: 0.0
wind_reference_y: 0.0
wind_reference_z: 0.0
wind_time_offset: 0.0
wind_enforce_bx_relation: 0
```

Keep the HDF5 beside that `config.yaml`, or change only `wind_file` to its
path relative to the config. The zeros are not placeholders:

- OMNI time is already referenced to bow-shock arrival, so the reference is
  native `(0,0,0)`;
- `/T=0` is exactly the requested simulation start, so the code-time offset
  is exactly zero;
- the measured OMNI `Bx` is retained rather than replaced by a fitted planar
  relation.

No student should calculate `4500/63.71`, guess a monitor position, or shift
the CSV in a spreadsheet. The model still applies its documented
position-dependent propagation from the bow-shock reference to each boundary
or grid point (section 6).

### 5.5 Use another event

Only change the UTC start, stop/duration, and output name:

```bash
python3 scripts/solar_wind/fetch_nasa_omni_event.py \
  --simulation-start START_UTC \
  --simulation-stop STOP_UTC \
  --output run_inputs/EVENT_NAME/EVENT_NAME_omni_1min.h5
```

Useful policy options are `--padding-minutes`, `--max-gap-minutes`, and
`--overwrite`; inspect all choices with `--help`. A new event can legitimately
fail if NASA has a longer outage than the declared policy or if its padding is
too short for the observed speed and 200-RE domain. That is a data-quality
decision, not a reason to bypass the checks. A network-free tool check is
available as:

```bash
python3 scripts/solar_wind/fetch_nasa_omni_event.py --self-test
```

If a collaborator already downloaded the exact headerless HAPI CSV, pass it
with `--raw-csv FILE`; the same validation and conversion path is used without
a second download. Raw ACE or DSCOVR spacecraft data are intentionally outside
this lazy event path because they require a separately justified propagation
analysis. Use them only when L1 propagation itself is part of the research.

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

For a NASA event, first open the generated `.qa.png` and
`.receipt.json`. Confirm that the blue-shaded simulation interval has the
expected IMF, velocity, density, and temperature, and that the receipt reports
an acceptable interpolation count. The one-command tool has already checked
the HDF5 schema, 60-second cadence, all primary values, sign mapping, coverage,
and source hashes.

An independent command-line inspection before an expensive run remains good
practice:

```bash
h5dump -n upstream_input.h5
h5dump -d /T -d /D -d /Vx -d /Bx -d /By -d /Bz upstream_input.h5
sha256sum upstream_input.h5
```

On macOS, `shasum -a 256 upstream_input.h5` is equivalent. Verify:

- one-dimensional equal-length arrays;
- finite values and strictly increasing time;
- positive density and thermal pressure or temperature;
- expected native `Vx` sign and no zero `Vx`;
- expected IMF signs after coordinate conversion;
- no source fill values or unexplained discontinuities;
- gap masks and interpolation policy agree with the generated QA/receipt;
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
  they were raw L1 observations. The checked-in NASA tool already prevents
  this; do not alter its generated reference or offset.
- **Manual NASA table assembly:** downloading eight separate columns, joining
  them, or shifting them in a spreadsheet bypasses the checked fill-value,
  frame, gap, and provenance path. Use `fetch_nasa_omni_event.py` once.
- **Hidden endpoint clamping:** a run appears to work but spends part of its
  interval at a frozen first or last sample.
- **Wrong interpolation:** `step` is appropriate for a designed discontinuous
  schedule, while observed event data normally require `linear`.
- **Partial coefficient override:** `wind_by_coefficient` and
  `wind_bz_coefficient` must be provided together.
- **Untracked preprocessing:** an event HDF5 without source, frame, epoch, gap,
  and hash metadata is not reproducible.
