#!/usr/bin/env bash
set -euo pipefail

source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-$source_dir/build/local}"
platform_cmake_args=()

if [[ "$(uname -s)" == "Darwin" ]]; then
  command -v brew >/dev/null || {
    echo "Homebrew is required; see $source_dir/docs/LOCAL_SETUP.md" >&2
    exit 2
  }
  mpi_root="$(brew --prefix open-mpi)"
  if brew list --versions hdf5-mpi >/dev/null 2>&1; then
    hdf5_root="$(brew --prefix hdf5-mpi)"
  else
    hdf5_root="$(brew --prefix hdf5)"
    echo "Using serial Homebrew HDF5 for the two-rank local smoke run."
    echo "Sugon production still requires parallel HDF5."
  fi
  c_compiler="$mpi_root/bin/mpicc"
else
  c_compiler="${MPICC:-$(command -v mpicc)}"
  hdf5_root="${HDF5_ROOT:-}"
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  # Homebrew OpenMPI exposes one local slot on some laptops. The smoke and
  # two-rank regression tests intentionally run two lightweight processes.
  platform_cmake_args+=("-DMPIEXEC_PREFLAGS=--map-by;slot:OVERSUBSCRIBE")
fi

cmake_args=(
  -S "$source_dir"
  -B "$build_dir"
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_C_COMPILER="$c_compiler"
  -DCOORD_TYPE=3
  -DPROBLEM=22
  -DHDF5_PREFER_PARALLEL=TRUE
  -DGAMERA_BOW_SHOCK_INNER_RADIUS_RE=2.5
  -DGAMERA_BOW_SHOCK_OUTER_RADIUS_RE=200.0
  -DGAMERA_EARTH_RADIAL_MAP_VERSION=4
  -DGAMERA_BOW_SHOCK_RADIAL_STRETCH=3.0
  -DGAMERA_STRICT_INNER_MAGNETIC_WALL=ON
  -DGAMERA_STRICT_INNER_GAS_WALL=ON
  -DGAMERA_YINYANG_MFE_INTERFACE=ON
  -DGAMERA_YINYANG_HDIV_RECONCILE=OFF
  -DGAMERA_YINYANG_HDIV_OPTIMIZED=OFF
  -DGAMERA_YINYANG_HDIV_DISTRIBUTED=OFF
  -DGAMERA_YINYANG_HDIV_PROFILE=OFF
  -DGAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY=OFF
  -DGAMERA_YINYANG_SPARSE_OVERSET=ON
  -DGAMERA_YINYANG_SPARSE_OVERSET_PROFILE=OFF
  -DGAMERA_BUILD_INNER_MAGNETOSPHERE_HPLUS=OFF
  -DGAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC=OFF
  -DGAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK=OFF
  -DGAMERA_BUILD_NONORTHOGONAL_TESTS=ON
)
if [[ -n "$hdf5_root" ]]; then
  cmake_args+=("-DHDF5_ROOT=$hdf5_root")
fi
cmake_args+=("${platform_cmake_args[@]}")

cmake "${cmake_args[@]}"
cmake --build "$build_dir" --parallel
ctest --test-dir "$build_dir" --output-on-failure

echo "Built and tested: $build_dir/run_mhd"
