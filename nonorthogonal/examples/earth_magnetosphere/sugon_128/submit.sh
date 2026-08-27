#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || "$1" != /* ]]; then
  echo "Usage: $0 /absolute/path/for/this/run" >&2
  exit 2
fi
run_root="$1"
example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

build_job="$(
  sbatch --parsable \
    --export="ALL,GAMERA_RUN_ROOT=$run_root" \
    "$example_dir/build.slurm"
)"
run_job="$(
  sbatch --parsable \
    --dependency="afterok:$build_job" \
    --export="ALL,GAMERA_RUN_ROOT=$run_root" \
    "$example_dir/run.slurm"
)"
diagnostics_job="$(
  sbatch --parsable \
    --dependency="afterok:$run_job" \
    --export="ALL,GAMERA_RUN_ROOT=$run_root" \
    "$example_dir/diagnostics.slurm"
)"

printf 'build_job=%s\nrun_job=%s\ndiagnostics_job=%s\n' \
  "$build_job" "$run_job" "$diagnostics_job"
