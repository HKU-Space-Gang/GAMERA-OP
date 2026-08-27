#!/bin/bash -l
#SBATCH --job-name=gamera-earth-build
#SBATCH --partition=hfacnormal01
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT
#SBATCH --qos=REPLACE_WITH_YOUR_QOS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue

set -euo pipefail

example_dir="$SLURM_SUBMIT_DIR"
source_dir="$(cd "$example_dir/../../.." && pwd)"
run_root="${GAMERA_RUN_ROOT:?submit with GAMERA_RUN_ROOT=/absolute/run/root}"
build_dir="$run_root/build"
cmake_bin="${CMAKE_BIN:-/public/home/binzheng/cmake-4.3.1-linux-x86_64/bin/cmake}"

module --force purge
module load compiler/gcc/12.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/hdf5/intelmpi/1.12.0
export I_MPI_CC="$(command -v gcc)"
export HDF5_ROOT="$HDF5_DIR"

test -d "$source_dir/src/nonorthogonal"
test -f "$example_dir/config.yaml"
test -f "$example_dir/../common/upstream_zero_then_minus5.h5"
test ! -e "$build_dir"
mkdir -p "$run_root"

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$(command -v mpicc)" \
  -DCOORD_TYPE=3 \
  -DPROBLEM=22 \
  -DHDF5_ROOT:PATH="$HDF5_ROOT" \
  -DHDF5_PREFER_PARALLEL=TRUE \
  -DGAMERA_BOW_SHOCK_INNER_RADIUS_RE=2.5 \
  -DGAMERA_BOW_SHOCK_OUTER_RADIUS_RE=200.0 \
  -DGAMERA_EARTH_RADIAL_MAP_VERSION=4 \
  -DGAMERA_BOW_SHOCK_RADIAL_STRETCH=3.0 \
  -DGAMERA_STRICT_KAIJU_MAGNETIC_WALL=ON \
  -DGAMERA_STRICT_KAIJU_GAS_WALL=ON \
  -DGAMERA_YINYANG_MFE_INTERFACE=ON \
  -DGAMERA_YINYANG_HDIV_RECONCILE=OFF \
  -DGAMERA_YINYANG_HDIV_OPTIMIZED=OFF \
  -DGAMERA_YINYANG_HDIV_DISTRIBUTED=OFF \
  -DGAMERA_YINYANG_HDIV_PROFILE=OFF \
  -DGAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY=OFF \
  -DGAMERA_YINYANG_SPARSE_OVERSET=ON \
  -DGAMERA_YINYANG_SPARSE_OVERSET_PROFILE=OFF \
  -DGAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS=OFF \
  -DGAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC=OFF \
  -DGAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK=OFF \
  -DGAMERA_BUILD_NONORTHOGONAL_TESTS=ON
"$cmake_bin" --build "$build_dir" --parallel 64

ctest --test-dir "$build_dir" --output-on-failure \
  -R '(nonorthogonal_(radial_map$|mi$|electron_precipitation$|mi_dpb$|mi_snapshot$|mi_collective_2rank$|absolute_cadence$|sha256$|mi_restart_2rank$)|yinyang_sparse_plan)'

test -x "$build_dir/run_mhd"
cp "$build_dir/run_mhd" "$run_root/run_mhd"
cp "$example_dir/config.yaml" "$run_root/config.yaml"
cp "$example_dir/../common/upstream_zero_then_minus5.h5" "$run_root/"
chmod 0750 "$run_root/run_mhd"

sha256sum "$run_root/run_mhd" "$run_root/config.yaml" \
  "$run_root/upstream_zero_then_minus5.h5" > "$run_root/build_inputs.sha256"
{
  "$cmake_bin" --version | head -1
  "$(command -v mpicc)" -show
  "$HDF5_ROOT/bin/h5pcc" -showconfig | grep -E 'Parallel HDF5:.*yes'
  module list 2>&1
} > "$run_root/build_environment.txt"
printf 'BUILD_COMPLETE binary=%s time=%s\n' \
  "$run_root/run_mhd" "$(date -Is)" | tee "$run_root/build_complete.txt"
