#!/bin/bash -l
#SBATCH --job-name=gamera-earth-128
#SBATCH --partition=hfacnormal01
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT
#SBATCH --qos=REPLACE_WITH_YOUR_QOS
#SBATCH --nodes=4
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=node380

set -euo pipefail

run_root="${GAMERA_RUN_ROOT:?submit with GAMERA_RUN_ROOT=/absolute/run/root}"
result_dir="$run_root/result"

module --force purge
module load compiler/gcc/12.3.0
module load mpi/intelmpi/2021.3.0
module load mathlib/hdf5/intelmpi/1.12.0

test -x "$run_root/run_mhd"
test -f "$run_root/config.yaml"
test -f "$run_root/upstream_zero_then_minus5.h5"
test ! -e "$result_dir"
mkdir -m 750 "$result_dir"
cp "$run_root/config.yaml" "$run_root/upstream_zero_then_minus5.h5" \
  "$result_dir/"
cd "$result_dir"

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export I_MPI_PIN=0
export HDF5_USE_FILE_LOCKING=FALSE

{
  echo "job_id=$SLURM_JOB_ID"
  echo "model=frozen_nonorthogonal_earth_production_v1"
  echo "grid=128x48x144x2"
  echo "mpi_ranks=256"
  echo "proc_dims_per_patch=4x4x8"
  echo "radial_map=v4_2p5_to_14_RE_uniform_0p5_RE_then_C1_exponential_to_200_RE"
  echo "mi_grid=384x43_per_hemisphere"
  echo "mi_solver=cached_conductance_tensor_BiCGStab_SSOR"
  echo "yinyang_sparse_overset=ON_profile=OFF"
  echo "precipitation=Fedder_Kaiju_beta_0.8724646"
  echo "euv=Kaiju_LOMPE_F10.7_150_twilight_tail_to_SZA_120"
  echo "sunward_direction=(-1,0,0)_Earth_to_Sun"
  echo "mfe=ON_hdiv=OFF_ring_current=OFF"
  sha256sum "$run_root/run_mhd" config.yaml upstream_zero_then_minus5.h5
  date -Is
} > run_receipt.txt

srun --mpi=pmi2 --cpu-bind=cores --kill-on-bad-exit=1 \
  "$run_root/run_mhd" config.yaml
printf 'RUN_COMPLETE rc=0 time=%s\n' "$(date -Is)" | tee run_complete.txt
