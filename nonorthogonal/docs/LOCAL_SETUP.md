# Local setup and first run

These instructions assume no prior GAMERA-OP experience. Commands are run in a
terminal. Do not type the leading `$` shown in some tutorials; the commands
below contain no prompt character.

## 1. Download the repository

```bash
git clone https://github.com/HKU-Space-Gang/GAMERA-OP.git
cd GAMERA-OP/nonorthogonal
```

Until the integration pull request is merged, replace the second command by
checking out the PR branch specified in that pull request.

## 2. Install compilers and libraries

The code needs a C compiler, CMake, MPI, OpenMP, HDF5, Python 3, NumPy, and
h5py. Parallel HDF5 is mandatory for production, but the small local smoke
test may use serial HDF5 because its compact files have a safe serial path.

### macOS with Homebrew

First install Homebrew from <https://brew.sh/> if `brew --version` fails. Then:

```bash
brew install cmake open-mpi hdf5 libomp python
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy h5py
```

On later terminal sessions, return to this directory and reactivate the
environment with `source .venv/bin/activate`.

### Ubuntu or Debian

```bash
sudo apt update
sudo apt install build-essential cmake openmpi-bin libopenmpi-dev \
  libhdf5-openmpi-dev libomp-dev python3 python3-venv
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy h5py
```

If CMake selects serial HDF5 on Ubuntu, set the parallel installation before
building:

```bash
export HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi
```

## 3. Verify the tools

```bash
cmake --version
mpicc --version
mpiexec --version
python3 -c "import h5py, numpy; print(h5py.__version__, numpy.__version__)"
```

All four commands must finish without an error.

## 4. Build and run the test suite

```bash
./scripts/build_local.sh
```

This creates `build/local/run_mhd` and runs the production regression tests.
The build command fixes the same numerical switches used for the accepted
Sugon run. It does not alter the standard solver in the repository root.

If you need a separate build directory:

```bash
./scripts/build_local.sh build/my_machine
```

## 5. Run the complete two-rank smoke case

```bash
./scripts/run_local_smoke.sh
```

This copies the input files into `run/local_smoke`, launches exactly two MPI
ranks, and checks the latest rank-local MHD restart and schema-3 MI HDF5 fields
for finite values and frozen-model metadata. Success ends with
`SMOKE_VALIDATION_PASS`. The laptop example uses restart output because compact
per-patch MHD analysis intentionally requires parallel HDF5; the Sugon
production example enables that path.

The smoke case starts at the Bz=-5 nT portion of the supplied wind history so
that the solar-wind reader and complete MI/precipitation/conductance path are
exercised immediately. It is an installation test, not a converged
magnetosphere simulation.

## 6. Common failures

### CMake found serial HDF5

For a production build, use `h5pcc -showconfig` and look for
`Parallel HDF5: yes`. Homebrew treats `hdf5` and `hdf5-mpi` as conflicting
variants; do not unlink a shared local installation just for the smoke test.
The supplied local script accepts serial HDF5, while the Sugon build script
explicitly requires the parallel module.

### The number of MPI ranks is wrong

The required count is twice the product of the three `proc_dims_*` values.
The local example is `2 * 1 * 1 * 1 = 2`.

### A run directory already exists

The script deliberately refuses to overwrite results. Choose a new directory:

```bash
./scripts/run_local_smoke.sh run/local_smoke_02
```

### macOS reports an OpenMP library error

Run `brew reinstall libomp`, then rebuild in a new build directory.

### MPI cannot start processes on a laptop

Confirm `mpiexec -n 2 hostname` works. Corporate security software can block
local MPI sockets; allow the MPI launcher and retry.
