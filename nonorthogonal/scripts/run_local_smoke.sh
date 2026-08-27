#!/usr/bin/env bash
set -euo pipefail

source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-$source_dir/build/local}"
run_dir="${1:-$source_dir/run/local_smoke}"
example_dir="$source_dir/examples/earth_magnetosphere"

test -x "$build_dir/run_mhd" || {
  echo "Missing $build_dir/run_mhd; run scripts/build_local.sh first." >&2
  exit 2
}
test ! -e "$run_dir" || {
  echo "Refusing to overwrite existing run directory: $run_dir" >&2
  exit 2
}

mkdir -p "$run_dir"
run_dir="$(cd "$run_dir" && pwd)"
cp "$example_dir/local_smoke/config.yaml" "$run_dir/"
cp "$example_dir/common/upstream_zero_then_minus5.h5" "$run_dir/"
cd "$run_dir"

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export HDF5_USE_FILE_LOCKING=FALSE

mpi_extra=()
if [[ "$(uname -s)" == "Darwin" ]]; then
  mpi_extra+=(--map-by slot:OVERSUBSCRIBE)
fi
mpiexec "${mpi_extra[@]}" -n 2 "$build_dir/run_mhd" config.yaml \
  2>&1 | tee run.log
"${PYTHON_BIN:-python3}" "$source_dir/scripts/validate_smoke.py" "$run_dir"
echo "Local smoke run passed: $run_dir"
