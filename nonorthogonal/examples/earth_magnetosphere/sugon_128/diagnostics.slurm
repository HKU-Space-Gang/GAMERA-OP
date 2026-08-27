#!/bin/bash -l
#SBATCH --job-name=gamera-earth-diag
#SBATCH --partition=hfacnormal01
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT
#SBATCH --qos=REPLACE_WITH_YOUR_QOS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --no-requeue
#SBATCH --exclude=node380

set -euo pipefail

example_dir="$SLURM_SUBMIT_DIR"
source_dir="$(cd "$example_dir/../../.." && pwd)"
diagnostics_dir="$source_dir/scripts/diagnostics"
run_root="${GAMERA_RUN_ROOT:?submit with GAMERA_RUN_ROOT=/absolute/run/root}"
result_dir="$run_root/result"
output_dir="$run_root/diagnostics"
python_bin="${PYTHON_BIN:-/public/software/apps/anaconda3/2023.09/bin/python3}"

test -f "$result_dir/run_complete.txt"
test -f "$result_dir/analysis_grid_p0.h5"
test -f "$diagnostics_dir/render_bzminus5_two_panel.py"
test ! -e "$output_dir"
mkdir -m 750 -p "$output_dir/logs" "$output_dir/tmp/matplotlib"

export MPLBACKEND=Agg
export MPLCONFIGDIR="$output_dir/tmp/matplotlib"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

(cd "$diagnostics_dir" && sha256sum -c DIAGNOSTICS_SHA256SUMS)

"$python_bin" "$diagnostics_dir/extract_radial_edges.py" \
  "$result_dir/analysis_grid_p0.h5" "$output_dir/radial_edges.json"

"$python_bin" "$diagnostics_dir/render_bzminus5_two_panel.py" \
  --inputs "$result_dir" \
  --output "$output_dir/standard" \
  --radial-edges "$output_dir/radial_edges.json" \
  --workers 8 \
  --start-seconds 0 \
  --end-seconds 18000 \
  --cadence-seconds 120 \
  --product-prefix gamera_earth_128 \
  > "$output_dir/logs/standard_render.log" 2>&1

final_mi="$(find "$result_dir" -maxdepth 1 -type f \
  -name 'mi_ionosphere_*.h5' -print | sort | tail -n 1)"
test -n "$final_mi"
"$python_bin" "$diagnostics_dir/render_mi_full_snapshot.py" \
  --input "$final_mi" \
  --output "$output_dir/mi_full_final.png" \
  > "$output_dir/logs/mi_snapshot.log" 2>&1

test "$(find "$output_dir/standard/frames" -maxdepth 1 -type f \
  -name 'frame_*.png' | wc -l)" -eq 151
test -s "$output_dir/standard/products/gamera_earth_128_two_panel_final.png"
test -s "$output_dir/standard/products/gamera_earth_128_mi_diagnostics.png"
test -s "$output_dir/mi_full_final.png"
tar -C "$output_dir/standard" -cf "$output_dir/standard_two_panel_frames.tar" frames

ffmpeg_bin="${FFMPEG_BIN:-$(command -v ffmpeg || true)}"
if [[ -n "$ffmpeg_bin" ]]; then
  "$ffmpeg_bin" -y -framerate 10 \
    -i "$output_dir/standard/frames/frame_%05d.png" \
    -c:v libx264 -pix_fmt yuv420p -movflags +faststart \
    "$output_dir/gamera_earth_128_two_panel.mp4" \
    > "$output_dir/logs/ffmpeg.log" 2>&1
fi

printf 'DIAGNOSTICS_COMPLETE frames=151 final_mi=%s time=%s\n' \
  "$final_mi" "$(date -Is)" | tee "$output_dir/diagnostics_complete.txt"
