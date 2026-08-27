# Production provenance

This release is numerically equivalent to the accepted source snapshot for the
Sugon five-hour Bz=-5 nT production run. Equations, constants, control flow and
production defaults are preserved. Public identifiers, comments, log labels
and test labels were normalized for this standalone GAMERA-OP release;
packaging also adds documentation, examples, scripts and a CMake fail-fast
guard while excluding unfinished research modules.

## Accepted run

- Slurm production job: `27542153`, `COMPLETED/0:0`
- elapsed time: `07:00:37`
- grid: `128x48x144x2`
- MPI ranks: `256`, decomposition `4x4x8` per patch
- physical sequence: zero IMF from 0 to 1 hour, then Bz=-5 nT through 5 hours
- MI grid: `384x43` per hemisphere
- MI cadence: 10 physical seconds
- compact MHD and schema-3 MI output cadence: 120 physical seconds
- restart cadence: 3600 physical seconds
- final MI time: 18,000 s
- final CPCP north/south: 112.748 / 116.4 kV
- final `max_divB_total`: 1.2725e-09
- complete `checkpoint_000005`: 256 rank files

The historical accepted-run executable SHA-256 was:

```text
86eee1e7fd7b2cbf509abd28eabde16afa88b3a6fb358c059c15181c3347b25b
```

This hash identifies the executable used for job `27542153`. A newly compiled
public-release executable has a different byte hash because public symbols and
log strings were renamed; regression tests verify the unchanged numerical
behavior. Current release source hashes are recorded in `SOURCE_SHA256SUMS`.

## Performance gates

The frozen build uses sparse Yin-Yang overset communication with profiling
disabled. A strict 120-second, 256-rank gate measured:

- 1032 steps in 150.44 solver-wall seconds;
- 0.1458 solver-wall seconds per step;
- 12.14 million cell-updates per second;
- simulation-seconds / wall-second = 0.7977.

MI uses a cached conductance tensor with BiCGStab/SSOR. In the strict 384x43
two-hemisphere gate, three complete MI updates took 0.439, 0.545 and 0.460 s;
the steady MI overhead was below one percent of the MHD work between coupling
updates.

## Important endpoint detail

The accepted run wrote 152 raw p0 and 152 raw p1 files. The last pair had
physical times 18,000.0 s and 18,000.000000000004 s. This is a non-physical
floating-point endpoint duplicate. It was retained as evidence. There are 151
unique 120-second physical times and exactly 151 MI files. Consumers must use
physical-time matching rather than raw filename counts.

## Frozen compile choices

- `COORD_TYPE=3`, `PROBLEM=22`
- inner radius 2.5 RE, outer radius 200 RE, radial map version 4
- strict production gas and magnetic inner walls ON
- MFE Yin-Yang interface ON
- sparse Yin-Yang overset ON, sparse profiling OFF
- H(div) reconciliation OFF
- inner-magnetosphere H+, ring current and angular-momentum correction OFF

The example configuration and input file are copied from the accepted run.
Their repository hashes are recorded in `SOURCE_SHA256SUMS` at release time.
