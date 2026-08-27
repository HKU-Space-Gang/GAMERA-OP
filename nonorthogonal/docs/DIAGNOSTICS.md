# Standard diagnostics and movie

These are the frozen Python renderers used for the accepted 128-resolution
production result. They are kept together in `scripts/diagnostics/`; do not
replace them with a case-specific plotting script when comparing production
runs.

## Included scripts

- `extract_radial_edges.py`: reads the exact radial vertices from
  `analysis_grid_p0.h5`.
- `postprocess_compact_base.py`: common compact-HDF discovery, restart
  deduplication, Yin-Yang interpolation and finite-field validation.
- `render_bzminus5_two_panel.py`: standard XY residual-Bz and XZ thermal-
  pressure figure, with embedded North/South FAC and potential panels. It
  creates one PNG per time, a final PNG, MHD and MI time-series PNGs, a CSV and
  a JSON manifest.
- `render_mi_full_snapshot.py`: North/South polar snapshots of potential,
  upward FAC, total Pedersen and Hall conductance, electron energy flux FE,
  and electron number flux FN.

`DIAGNOSTICS_SHA256SUMS` records the accepted copies. From the diagnostics
directory, `sha256sum -c DIAGNOSTICS_SHA256SUMS` must pass.

## Python requirements

Activate the environment described in `LOCAL_SETUP.md`, then install:

```bash
python3 -m pip install numpy h5py matplotlib
```

Install `ffmpeg` only if an MP4 is required. On macOS use
`brew install ffmpeg`; on Ubuntu use `sudo apt install ffmpeg`.

Set a non-interactive backend on a cluster:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/absolute/writable/path/matplotlib-cache
```

## Expected production files

In the commands below, `RESULT` is a completed run directory containing:

```text
analysis_grid_p0.h5       analysis_grid_p1.h5
analysis_p0_*.h5          analysis_p1_*.h5
mi_ionosphere_*.h5
```

Compact MHD output requires parallel HDF5. The two-rank laptop smoke example
uses rank-local restart files instead and therefore is not an input to the
standard MHD renderer.

## Render all standard frames and diagnostics

From `GAMERA-OP/nonorthogonal`:

```bash
RESULT=/absolute/path/to/completed/result
OUT=/absolute/path/to/fresh/diagnostics
PYTHON=python3

mkdir -p "$OUT"
$PYTHON scripts/diagnostics/extract_radial_edges.py \
  "$RESULT/analysis_grid_p0.h5" "$OUT/radial_edges.json"

$PYTHON scripts/diagnostics/render_bzminus5_two_panel.py \
  --inputs "$RESULT" \
  --output "$OUT/standard" \
  --radial-edges "$OUT/radial_edges.json" \
  --workers 8 \
  --start-seconds 0 \
  --end-seconds 18000 \
  --cadence-seconds 120 \
  --product-prefix gamera_earth_128
```

The output directory must be fresh; the renderer refuses to overwrite an
existing `frames/` or `products/` directory. For the supplied five-hour case,
151 frames are expected. To preview only one physical frame, add
`--frame-index 150` and use a different fresh output directory.

Important products are:

```text
standard/frames/frame_00000.png ... frame_00150.png
standard/products/gamera_earth_128_two_panel_final.png
standard/products/gamera_earth_128_mhd_x0_y20_diagnostics.png
standard/products/gamera_earth_128_mi_diagnostics.png
standard/products/gamera_earth_128_diagnostics.csv
standard/products/gamera_earth_128_two_panel_manifest.json
```

The two-panel renderer is frozen for the accepted 128 radial-cell production
grid. Its time range and cadence are command-line inputs; changing the radial
cell count requires a separately reviewed renderer change.

## Render a full MI snapshot

For the accepted run the final MI file is `mi_ionosphere_001801.h5`. For a
different run, select the MI file whose root `time_seconds` attribute matches
the desired physical time.

```bash
$PYTHON scripts/diagnostics/render_mi_full_snapshot.py \
  --input "$RESULT/mi_ionosphere_001801.h5" \
  --output "$OUT/mi_full_t18000.png"
```

This script reads both hemispheres from the same schema-3 HDF5 file. It does
not mirror the South longitude axis. Both panels use 12 MLT at top, 06 at
right, 00 at bottom and 18 at left.

## Encode the movie

```bash
ffmpeg -y -framerate 10 \
  -i "$OUT/standard/frames/frame_%05d.png" \
  -c:v libx264 -pix_fmt yuv420p -movflags +faststart \
  "$OUT/gamera_earth_128_two_panel.mp4"

ffprobe -v error -select_streams v:0 \
  -show_entries stream=codec_name,pix_fmt,nb_frames,r_frame_rate \
  -of default=noprint_wrappers=1 \
  "$OUT/gamera_earth_128_two_panel.mp4"
```

The expected metadata are H.264, `yuv420p`, 151 frames and 10 frames/s.

## Restart and endpoint handling

`postprocess_compact_base.py` groups p0 and p1 by rounded physical seconds,
checks patch-time agreement, and selects one pair per requested cadence. This
is intentional: a run that lands on the floating-point endpoint may contain
two raw MHD pairs at effectively 18,000 s. Keep both raw files as provenance;
the renderer must produce one physical frame.

Multiple directories may be passed after `--inputs` to join restart segments.
For duplicate physical times, the later input directory has priority.
