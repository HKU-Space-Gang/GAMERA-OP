#!/usr/bin/env python3
"""Small integrity gate for the packaged Earth magnetosphere smoke run."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

try:
    import h5py
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "Validation needs numpy and h5py. Install them with: "
        "python3 -m pip install numpy h5py"
    ) from exc


def assert_finite(path: Path) -> None:
    with h5py.File(path, "r") as handle:
        def visit(name: str, obj: h5py.Dataset) -> None:
            if isinstance(obj, h5py.Dataset) and np.issubdtype(obj.dtype, np.number):
                values = obj[...]
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"non-finite dataset {name} in {path}")

        handle.visititems(visit)


def assert_mi_science_finite(path: Path) -> None:
    science_fields = (
        "DPB_boundary_MLAT_deg",
        "DPB_mask",
        "F_E_erg_cm2_s",
        "F_N_cm2_s",
        "Sigma_H_EUV_S",
        "Sigma_H_S",
        "Sigma_P_EUV_S",
        "Sigma_P_S",
        "diffuse_precipitation_selector",
        "fac_parallel_A_m2",
        "potential_V",
    )
    with h5py.File(path, "r") as handle:
        for name in ("longitude_rad", "colatitude_rad"):
            if not np.all(np.isfinite(handle[name][...])):
                raise ValueError(f"non-finite MI coordinate {name} in {path}")
        for hemisphere in ("north", "south"):
            for name in science_fields:
                values = handle[f"{hemisphere}/{name}"][...]
                if not np.all(np.isfinite(values)):
                    raise ValueError(
                        f"non-finite MI science field {hemisphere}/{name} in {path}"
                    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    run_dir = args.run_directory.resolve()

    p0 = sorted(run_dir.glob("analysis_p0_*.h5"))
    p1 = sorted(run_dir.glob("analysis_p1_*.h5"))
    output_kind = "compact-analysis"
    if not p0 and not p1:
        p0 = sorted(run_dir.rglob("restart_p0_*.h5"))
        p1 = sorted(run_dir.rglob("restart_p1_*.h5"))
        output_kind = "rank-local-restart"
    mi = sorted(run_dir.glob("mi_ionosphere_*.h5"))
    if not p0 or len(p0) != len(p1) or not mi:
        raise SystemExit(
            f"incomplete output: p0={len(p0)}, p1={len(p1)}, mi={len(mi)}"
        )

    for path in (p0[-1], p1[-1]):
        assert_finite(path)
    assert_mi_science_finite(mi[-1])

    with h5py.File(mi[-1], "r") as handle:
        expected = {
            "schema_version": 3,
            "hemisphere_longitude_axes_identical": 1,
        }
        for name, value in expected.items():
            if int(handle.attrs.get(name, -1)) != value:
                raise ValueError(f"{name} is not {value} in {mi[-1]}")
        beta = float(handle.attrs.get("electron_precipitation_beta", np.nan))
        sunward = np.array(
            [
                handle.attrs.get("sunward_direction_x", np.nan),
                handle.attrs.get("sunward_direction_y", np.nan),
                handle.attrs.get("sunward_direction_z", np.nan),
            ],
            dtype=float,
        )
        if not np.isclose(beta, 0.8724646, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"unexpected precipitation beta: {beta}")
        if not np.allclose(
            sunward, [-1.0, 0.0, 0.0], rtol=0.0, atol=1.0e-12
        ):
            raise ValueError(f"unexpected Earth-to-Sun direction: {sunward}")

    print(
        "SMOKE_VALIDATION_PASS "
        f"mhd={output_kind} p0={len(p0)} p1={len(p1)} mi={len(mi)} "
        f"latest={mi[-1].name}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
