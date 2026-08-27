# Sugon production setup

The supplied example reproduces the accepted 128x48x144x2 Earth run on 256
MPI ranks. It uses four nodes with 64 ranks per node and one OpenMP thread per
rank.

## 1. Clone once on the login node

```bash
cd /public/home/YOUR_USER
git clone https://github.com/HKU-Space-Gang/GAMERA-OP.git
cd GAMERA-OP/nonorthogonal/examples/earth_magnetosphere/sugon_128
```

Use a Git tag or commit SHA for real production, not a moving branch.

## 2. Set Slurm account information

Edit `build.slurm` and `run.slurm`. Replace:

```text
REPLACE_WITH_YOUR_ACCOUNT
REPLACE_WITH_YOUR_QOS
```

If your allocation uses another partition or excludes different nodes, change
only those `#SBATCH` lines. The accepted software modules are:

```text
compiler/gcc/12.3.0
mpi/intelmpi/2021.3.0
mathlib/hdf5/intelmpi/1.12.0
```

The example defaults to CMake 4.3.1 at
`/public/home/binzheng/cmake-4.3.1-linux-x86_64/bin/cmake`. If that file is not
readable from your account, set an installed CMake explicitly when submitting:

```bash
export CMAKE_BIN="$(command -v cmake)"
```

## 3. Submit the dependency chain

Choose a fresh absolute run directory. The scripts refuse to overwrite it.

```bash
./submit.sh /public/home/YOUR_USER/runs/gamera_earth_bzminus5_128_v1
```

The command prints two IDs. The production run has an `afterok` dependency on
the build and regression gate, so a failed build cannot launch a large job.

Check status with:

```bash
squeue -j BUILD_JOB_ID,RUN_JOB_ID
sacct -j BUILD_JOB_ID,RUN_JOB_ID \
  --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS
```

Follow the production log with the output filename reported by Slurm, for
example:

```bash
tail -f slurm-RUN_JOB_ID.out
```

## 4. Acceptance checks

The build must finish `COMPLETED/0:0`, all selected regression tests must pass,
and `build_inputs.sha256` must exist in the run root.

The five-hour run must finish `COMPLETED/0:0` and contain:

- 151 schema-3 MI files at 120-second physical cadence;
- a complete `checkpoint_000005` with 256 rank files;
- finite final p0, p1, and MI fields;
- no ERROR, NaN, Inf, MPI exit, or failed MI convergence record;
- `run_receipt.txt` and `run_complete.txt`.

The accepted source can produce a second raw p0/p1 file at the final time that
differs from 18,000 s only by floating-point roundoff. Keep both raw files,
but count unique physical times and select the closest file for diagnostics.
There must still be exactly 151 unique 120-second physical times and 151 MI
files.

## 5. What may be changed for a science run

Make a copy of `config.yaml` in a new run directory and record every change.
Solar-wind history and runtime are normal science inputs. Do not change the
radial map, Yin-Yang interface, sparse-overset or MI solver build switches in
an A/B science comparison unless the numerical change itself is the research
question.

Never reuse a result directory. The supplied scripts make inputs immutable by
copying them into the result and recording SHA-256 hashes before launch.
