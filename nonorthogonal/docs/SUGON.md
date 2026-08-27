# Sugon production setup

The supplied example reproduces the accepted 128x48x144x2 Earth run on 256
MPI ranks. It uses four nodes with 64 ranks per node and one OpenMP thread per
rank.

## 1. Clone once on the login node

```bash
cd /public/home/YOUR_USER
git clone https://github.com/ijmhd/GAMERA-OP.git
cd GAMERA-OP/nonorthogonal/examples/earth_magnetosphere/sugon_128
```

Use a Git tag or commit SHA for real production, not a moving branch.

This checked-in example is the frozen 128 preset. For the exact differences
between the 32 smoke, 72 medium and 128 production grids, see
[Grid and MPI options](GRID_OPTIONS.md).

## 2. Set Slurm account information

Edit `build.slurm`, `run.slurm`, and `diagnostics.slurm`. Replace:

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

The command prints three IDs. The production run has an `afterok` dependency
on the build and regression gate, and the standard diagnostics job has an
`afterok` dependency on production. A failed stage therefore cannot launch a
later stage.

Check status with:

```bash
squeue -j BUILD_JOB_ID,RUN_JOB_ID,DIAGNOSTICS_JOB_ID
sacct -j BUILD_JOB_ID,RUN_JOB_ID,DIAGNOSTICS_JOB_ID \
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

The dependent `diagnostics.slurm` job uses the frozen scripts documented in
[Standard diagnostics and movie](DIAGNOSTICS.md). It must finish with
`DIAGNOSTICS_COMPLETE`, 151 non-empty frames, a final standard two-panel PNG,
MHD and MI time-series PNGs, a full North/South MI snapshot and a frames tar.
If `ffmpeg` is available, it also writes a 10-fps H.264 `yuv420p` MP4. Set
`FFMPEG_BIN=/absolute/path/to/ffmpeg` at submission if it is not on `PATH`.

The diagnostic products are placed under:

```text
RUN_ROOT/diagnostics/
```

For manual rendering, alternate time windows and the exact Python commands,
follow [Standard diagnostics and movie](DIAGNOSTICS.md). Dataset layouts and
coordinate conventions are specified in
[Yin-Yang output data model](OUTPUT_DATA_MODEL.md).

## 5. What may be changed for a science run

Make a copy of `config.yaml` in a new run directory and record every change.
Solar-wind history and runtime are normal science inputs. Use
[Solar-wind input](SOLAR_WIND.md) to construct and validate either an
idealized schedule or an observed event before submission. The supported NASA
event workflow is one Python command: it makes the raw download, HDF5, exact
YAML fragment, QA plot, and checksum receipt, so no manual time shift is
needed. Copy the resulting HDF5 beside the run `config.yaml` and use the
generated `.wind.yaml` settings. Do not change the
radial map, Yin-Yang interface, sparse-overset or MI solver build switches in
an A/B science comparison unless the numerical change itself is the research
question.

Never reuse a result directory. The supplied scripts make inputs immutable by
copying them into the result and recording SHA-256 hashes before launch.
