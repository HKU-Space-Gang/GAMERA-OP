#!/usr/bin/env python3
"""Restart-deduplicated diagnostics and movie frames for compact Yin-Yang output."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
from pathlib import Path
import shutil

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np


FIELDS = ("rho", "p", "vx", "vy", "vz", "bx", "by", "bz")
PROTON_MASS_KG = 1.67262192369e-27
EXPECTED_TIMES = np.arange(0.0, 21600.0 + 1.0, 120.0)
WORKER_CACHE: dict[str, object] = {}


def scalar(handle: h5py.File, name: str) -> float:
    return float(np.asarray(handle[name]).reshape(-1)[0])


def discover_mhd(inputs: list[Path]) -> list[dict[str, object]]:
    by_time: dict[int, dict[str, object]] = {}
    for priority, root in enumerate(inputs):
        for p0 in sorted(root.glob("analysis_p0_*.h5")):
            p1 = root / p0.name.replace("analysis_p0_", "analysis_p1_")
            if not p1.is_file():
                raise RuntimeError(f"missing Yang patch for {p0}")
            with h5py.File(p0, "r") as left, h5py.File(p1, "r") as right:
                time_seconds = scalar(left, "time_seconds")
                if abs(time_seconds - scalar(right, "time_seconds")) > 1.0e-6:
                    raise RuntimeError(f"patch time mismatch: {p0}, {p1}")
            key = int(round(time_seconds))
            if (key < int(round(EXPECTED_TIMES[0])) or
                    key > int(round(EXPECTED_TIMES[-1]))):
                continue
            by_time[key] = {
                "time_seconds": time_seconds,
                "p0": str(p0),
                "p1": str(p1),
                "source_priority": priority,
            }
    records = [by_time[key] for key in sorted(by_time)]
    if len(records) != len(EXPECTED_TIMES):
        raise RuntimeError(
            f"expected {len(EXPECTED_TIMES)} restart-deduplicated MHD frames, "
            f"found {len(records)}"
        )
    actual = np.asarray([record["time_seconds"] for record in records])
    if np.max(np.abs(actual - EXPECTED_TIMES)) > 1.0:
        raise RuntimeError(f"MHD cadence mismatch: max error {np.max(np.abs(actual - EXPECTED_TIMES))}")
    return records


def discover_mi(inputs: list[Path]) -> dict[int, str]:
    by_time: dict[int, str] = {}
    for root in inputs:
        for path in sorted(root.glob("mi_ionosphere_*.h5")):
            with h5py.File(path, "r") as handle:
                time_seconds = float(handle.attrs["time_seconds"])
            by_time[int(round(time_seconds))] = str(path)
    missing = [int(value) for value in EXPECTED_TIMES if int(round(value)) not in by_time]
    if missing:
        raise RuntimeError(f"missing M-I times: {missing[:10]}")
    return by_time


def logical_index(q, lower, upper, count):
    return (q - lower) * count / (upper - lower) - 0.5


def radial_index(radius, edges):
    centers = 0.5 * (edges[:-1] + edges[1:])
    return np.interp(radius, centers, np.arange(len(centers), dtype=float))


def patch_geometry(patch, x, y, z, lower, upper, shape, radial_edges):
    if patch == 0:
        lx, ly, lz = x, y, z
    else:
        lx, ly, lz = -x, z, y
    radius = np.sqrt(lx * lx + ly * ly + lz * lz)
    theta = np.arccos(np.clip(lz / radius, -1.0, 1.0))
    phi = np.arctan2(ly, lx)
    q = (
        radial_index(radius, radial_edges),
        logical_index(theta, lower[1], upper[1], shape[1]),
        logical_index(phi, lower[2], upper[2], shape[2]),
    )
    valid = np.ones(x.shape, dtype=bool)
    base, fraction = [], []
    for axis in range(3):
        valid &= (q[axis] >= -0.5) & (q[axis] <= shape[axis] - 0.5)
        clipped = np.clip(q[axis], 0.0, shape[axis] - 1.0)
        ibase = np.minimum(np.floor(clipped).astype(np.int64), shape[axis] - 2)
        base.append(ibase)
        fraction.append(clipped - ibase)
    margin = np.minimum.reduce(
        (q[1] + 0.5, shape[1] - 0.5 - q[1],
         q[2] + 0.5, shape[2] - 0.5 - q[2])
    )
    return valid, margin, tuple(base), tuple(fraction)


def interpolate(values, base, fraction):
    result = np.zeros(base[0].shape, dtype=np.float64)
    for di in (0, 1):
        wi = fraction[0] if di else 1.0 - fraction[0]
        for dj in (0, 1):
            wj = fraction[1] if dj else 1.0 - fraction[1]
            for dk in (0, 1):
                wk = fraction[2] if dk else 1.0 - fraction[2]
                result += wi * wj * wk * values[
                    base[0] + di, base[1] + dj, base[2] + dk
                ]
    return result


def interpolate_point(dataset, base, fraction) -> float:
    result = 0.0
    for di in (0, 1):
        wi = float(fraction[0][0]) if di else 1.0 - float(fraction[0][0])
        for dj in (0, 1):
            wj = float(fraction[1][0]) if dj else 1.0 - float(fraction[1][0])
            for dk in (0, 1):
                wk = float(fraction[2][0]) if dk else 1.0 - float(fraction[2][0])
                result += wi * wj * wk * float(dataset[
                    int(base[0][0]) + di,
                    int(base[1][0]) + dj,
                    int(base[2][0]) + dk,
                ])
    return result


def field_names(handle: h5py.File) -> list[str]:
    return [
        bytes(value).split(b"\0", 1)[0].decode("ascii").lower()
        for value in np.asarray(handle["field_names"])
    ]


def build_geometry(handle: h5py.File, radial_edges: np.ndarray) -> dict[str, object]:
    shape = tuple(int(value) for value in np.asarray(handle["global_cells"]))
    lower = np.asarray(handle["logical_lower"], dtype=float)
    upper = np.asarray(handle["logical_upper"], dtype=float)
    r = 0.5 * (radial_edges[:-1] + radial_edges[1:])
    nangle = 384

    phi_edges = np.linspace(-np.pi, np.pi, nangle + 1)
    phi = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    rr_xy, pp_xy = np.meshgrid(r, phi, indexing="ij")
    x_xy = rr_xy * np.cos(pp_xy)
    y_xy = rr_xy * np.sin(pp_xy)
    z_xy = np.zeros_like(x_xy)
    rr_xy_e, pp_xy_e = np.meshgrid(radial_edges, phi_edges, indexing="ij")

    alpha_edges = np.linspace(-np.pi, np.pi, nangle + 1)
    alpha = 0.5 * (alpha_edges[:-1] + alpha_edges[1:])
    rr_xz, aa_xz = np.meshgrid(r, alpha, indexing="ij")
    x_xz = rr_xz * np.cos(aa_xz)
    y_xz = np.zeros_like(x_xz)
    z_xz = rr_xz * np.sin(aa_xz)
    rr_xz_e, aa_xz_e = np.meshgrid(radial_edges, alpha_edges, indexing="ij")

    probe = (np.asarray([0.0]), np.asarray([20.0]), np.asarray([0.0]))
    return {
        "shape": shape,
        "lower": lower,
        "upper": upper,
        "x_xy": x_xy,
        "y_xy": y_xy,
        "z_xy": z_xy,
        "x_xy_e": rr_xy_e * np.cos(pp_xy_e),
        "y_xy_e": rr_xy_e * np.sin(pp_xy_e),
        "x_xz": x_xz,
        "y_xz": y_xz,
        "z_xz": z_xz,
        "x_xz_e": rr_xz_e * np.cos(aa_xz_e),
        "z_xz_e": rr_xz_e * np.sin(aa_xz_e),
        "plane_xy": [patch_geometry(patch, x_xy, y_xy, z_xy, lower, upper, shape, radial_edges)
                     for patch in (0, 1)],
        "plane_xz": [patch_geometry(patch, x_xz, y_xz, z_xz, lower, upper, shape, radial_edges)
                     for patch in (0, 1)],
        "probe": [patch_geometry(patch, *probe, lower, upper, shape, radial_edges)
                  for patch in (0, 1)],
    }


def choose_patch(geometries) -> np.ndarray:
    valid0, margin0 = geometries[0][:2]
    valid1, margin1 = geometries[1][:2]
    if np.any(~valid0 & ~valid1):
        raise RuntimeError("uncovered Yin-Yang target")
    return valid1 & (~valid0 | (margin1 > margin0))


def compose_field(handles, indices, geometries, name: str) -> np.ndarray:
    values = []
    for handle, index, geometry in zip(handles, indices, geometries):
        raw = np.asarray(handle["state"][index[name]], dtype=np.float64)
        values.append(interpolate(raw, geometry[2], geometry[3]))
    return np.where(choose_patch(geometries), values[1], values[0])


def probe_fields(handles, indices, geometries) -> dict[str, float]:
    selected = int(bool(choose_patch(geometries)[0]))
    geometry = geometries[selected]
    return {
        name: interpolate_point(
            handles[selected]["state"][indices[selected][name]],
            geometry[2], geometry[3],
        )
        for name in FIELDS
    }


def dipole_bz(x, y, z, moment_z):
    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    return moment_z * (3.0 * z * z / (r2 * r2 * r) - 1.0 / (r2 * r))


def cyclic(values, longitude):
    return (
        np.concatenate((longitude, [longitude[0] + 2.0 * np.pi])),
        np.concatenate((values, values[:, :1]), axis=1),
    )


def mi_arrays(path: str) -> tuple[dict[str, dict[str, object]], np.ndarray, np.ndarray]:
    result: dict[str, dict[str, object]] = {}
    with h5py.File(path, "r") as handle:
        longitude = np.asarray(handle["longitude_rad"], dtype=float)
        colatitude = np.asarray(handle["colatitude_rad"], dtype=float)
        for hemisphere in ("north", "south"):
            group = handle[hemisphere]
            fallback = -1 if hemisphere == "north" else 1
            multiplier = int(group.attrs.get("upward_fac_multiplier", fallback))
            fac = np.asarray(group["fac_parallel_A_m2"], dtype=float) * multiplier * 1.0e6
            potential = np.asarray(group["potential_V"], dtype=float) * 1.0e-3
            result[hemisphere] = {
                "fac": fac,
                "potential": potential,
                "cpcp_kv": float(np.ptp(potential)),
                "fac_min": float(np.min(fac)),
                "fac_max": float(np.max(fac)),
                "fac_rms": float(np.sqrt(np.mean(fac * fac))),
            }
    return result, longitude, colatitude


def polar_panel(axis, data, longitude, colatitude, hemisphere):
    lon, fac = cyclic(data["fac"], longitude)
    _, potential = cyclic(data["potential"], longitude)
    ll, cc = np.meshgrid(lon, colatitude)
    mesh = axis.pcolormesh(
        ll, cc, fac, shading="auto", cmap="RdBu_r",
        norm=colors.TwoSlopeNorm(vmin=-4.5, vcenter=0.0, vmax=4.5),
        rasterized=True,
    )
    levels = np.linspace(-70.0, 70.0, 15)
    contour = axis.contour(ll, cc, potential, levels=levels, colors="k", linewidths=0.45)
    axis.clabel(contour, contour.levels[::2], fontsize=6, fmt="%.0f")
    axis.set_theta_zero_location("N")
    axis.set_theta_direction(-1)
    axis.set_ylim(0.0, float(colatitude[-1]))
    axis.set_xticks(np.deg2rad([0, 90, 180, 270]))
    axis.set_xticklabels(["12", "06", "00", "18"], fontsize=8)
    axis.set_yticklabels([])
    axis.grid(alpha=0.35, linewidth=0.45)
    axis.set_title(
        f"{hemisphere.capitalize()} FAC / $\\Phi$\nCPCP={data['cpcp_kv']:.1f} kV",
        fontsize=9,
    )
    return mesh


def initialize_worker(radial_edges_list):
    WORKER_CACHE["radial_edges"] = np.asarray(radial_edges_list, dtype=float)


def render_frame(task: dict[str, object]) -> dict[str, object]:
    radial_edges = WORKER_CACHE["radial_edges"]
    paths = [task["p0"], task["p1"]]
    handles = [h5py.File(path, "r") for path in paths]
    try:
        names = [field_names(handle) for handle in handles]
        indices = [{name: values.index(name) for name in FIELDS} for values in names]
        if "geometry" not in WORKER_CACHE:
            WORKER_CACHE["geometry"] = build_geometry(handles[0], radial_edges)
        geometry = WORKER_CACHE["geometry"]
        normalization = np.asarray(handles[0]["normalization"], dtype=float)

        bz = compose_field(handles, indices, geometry["plane_xy"], "bz")
        pressure_code = compose_field(handles, indices, geometry["plane_xz"], "p")
        probe = probe_fields(handles, indices, geometry["probe"])
        moment_z = -3.1e-5 / normalization[6]
        residual_bz = (
            bz - dipole_bz(geometry["x_xy"], geometry["y_xy"], geometry["z_xy"], moment_z)
        ) * normalization[6] * 1.0e9
        pressure = pressure_code * normalization[5] * 1.0e9

        mi, longitude, colatitude = mi_arrays(task["mi"])
        figure = plt.figure(figsize=(15.5, 9.0), constrained_layout=True)
        grid = figure.add_gridspec(2, 3, width_ratios=(1.3, 1.3, 0.9))
        ax_xy = figure.add_subplot(grid[0, :2])
        ax_xz = figure.add_subplot(grid[1, :2])
        ax_n = figure.add_subplot(grid[0, 2], projection="polar")
        ax_s = figure.add_subplot(grid[1, 2], projection="polar")

        bmesh = ax_xy.pcolormesh(
            -geometry["x_xy_e"], geometry["y_xy_e"], residual_bz,
            shading="flat", cmap="RdBu_r",
            norm=colors.TwoSlopeNorm(vmin=-20.0, vcenter=0.0, vmax=20.0),
            rasterized=True,
        )
        pmesh = ax_xz.pcolormesh(
            -geometry["x_xz_e"], geometry["z_xz_e"], pressure,
            shading="flat", cmap="magma", norm=colors.LogNorm(vmin=0.01, vmax=10.0),
            rasterized=True,
        )
        for axis in (ax_xy, ax_xz):
            axis.add_patch(plt.Circle((0, 0), 3.0, color="0.25", zorder=5))
            axis.set_xlim(-145.0, 30.0)
            axis.set_ylim(-60.0, 60.0)
            axis.set_aspect("equal", adjustable="box")
            axis.grid(alpha=0.2, linewidth=0.4)
            axis.set_xlabel("SM X [$R_E$]")
        ax_xy.set_ylabel("SM Y [$R_E$]")
        ax_xz.set_ylabel("SM Z [$R_E$]")
        ax_xy.set_title(r"XY: $B_z-B_{z,\mathrm{dipole}}$")
        ax_xz.set_title("XZ: thermal pressure")
        figure.colorbar(bmesh, ax=ax_xy, pad=0.01, fraction=0.025).set_label("Residual $B_z$ [nT]")
        figure.colorbar(pmesh, ax=ax_xz, pad=0.01, fraction=0.025).set_label("Pressure [nPa]")

        fac_mesh = polar_panel(ax_n, mi["north"], longitude, colatitude, "north")
        polar_panel(ax_s, mi["south"], longitude, colatitude, "south")
        figure.colorbar(fac_mesh, ax=(ax_n, ax_s), pad=0.08, fraction=0.05).set_label(
            r"Upward FAC [$\mu$A m$^{-2}$]"
        )

        time_seconds = float(task["time_seconds"])
        if time_seconds < 3600.0:
            driving = "zero-IMF preconditioning"
        elif time_seconds < 10800.0:
            driving = r"northward IMF $B_z=+5$ nT"
        else:
            driving = r"southward IMF $B_z=-5$ nT"
        figure.suptitle(
            f"Yin–Yang Earth magnetosphere 256×96×288×2, 2048 MPI ranks — "
            f"t={time_seconds:.0f} s ({time_seconds / 3600.0:.2f} h), {driving}",
            fontsize=14,
        )
        output = Path(task["frame"])
        figure.savefig(output, dpi=120, facecolor="white")
        plt.close(figure)

        density_cm3 = probe["rho"] * normalization[4] / PROTON_MASS_KG / 1.0e6
        result: dict[str, object] = {
            "frame": int(task["frame_index"]),
            "time_seconds": time_seconds,
            "frame_path": str(output),
            "probe_density_cm3": density_cm3,
            "probe_pressure_npa": probe["p"] * normalization[5] * 1.0e9,
            "xy_residual_bz_min_nt": float(np.min(residual_bz)),
            "xy_residual_bz_max_nt": float(np.max(residual_bz)),
            "xz_pressure_min_npa": float(np.min(pressure)),
            "xz_pressure_max_npa": float(np.max(pressure)),
        }
        for name in ("vx", "vy", "vz"):
            result[f"probe_{name}_kms"] = probe[name] * normalization[2] / 1000.0
        for name in ("bx", "by", "bz"):
            result[f"probe_{name}_nt"] = probe[name] * normalization[6] * 1.0e9
        for hemisphere in ("north", "south"):
            result[f"cpcp_{hemisphere}_kv"] = mi[hemisphere]["cpcp_kv"]
            result[f"fac_{hemisphere}_min_uam2"] = mi[hemisphere]["fac_min"]
            result[f"fac_{hemisphere}_max_uam2"] = mi[hemisphere]["fac_max"]
            result[f"fac_{hemisphere}_rms_uam2"] = mi[hemisphere]["fac_rms"]
        return result
    finally:
        for handle in handles:
            handle.close()


def plot_diagnostics(records: list[dict[str, object]], products: Path) -> None:
    time_hours = np.asarray([record["time_seconds"] for record in records]) / 3600.0
    def values(name):
        return np.asarray([record[name] for record in records], dtype=float)

    figure, axes = plt.subplots(4, 1, figsize=(12.5, 12.0), sharex=True, constrained_layout=True)
    axes[0].plot(time_hours, values("probe_density_cm3"), label="density [cm$^{-3}$]")
    pressure_axis = axes[0].twinx()
    pressure_axis.plot(time_hours, values("probe_pressure_npa"), color="tab:red", label="pressure [nPa]")
    axes[0].set_ylabel("Density [cm$^{-3}$]")
    pressure_axis.set_ylabel("Pressure [nPa]")
    axes[1].plot(time_hours, values("probe_vx_kms"), label="$V_x$")
    axes[1].plot(time_hours, values("probe_vy_kms"), label="$V_y$")
    axes[1].plot(time_hours, values("probe_vz_kms"), label="$V_z$")
    axes[1].set_ylabel("Velocity [km/s]")
    axes[1].legend(ncol=3)
    axes[2].plot(time_hours, values("probe_bx_nt"), label="$B_x$")
    axes[2].plot(time_hours, values("probe_by_nt"), label="$B_y$")
    axes[2].plot(time_hours, values("probe_bz_nt"), label="$B_z$")
    axes[2].set_ylabel("Magnetic field [nT]")
    axes[2].legend(ncol=3)
    axes[3].plot(time_hours, values("xy_residual_bz_min_nt"), label="XY min residual $B_z$")
    axes[3].plot(time_hours, values("xy_residual_bz_max_nt"), label="XY max residual $B_z$")
    axes[3].set_ylabel("Residual $B_z$ [nT]")
    axes[3].set_xlabel("Physical time [h]")
    axes[3].legend(ncol=2)
    for axis in axes:
        axis.axvline(1.0, color="0.4", ls="--", lw=0.8)
        axis.axvline(3.0, color="0.4", ls="--", lw=0.8)
        axis.grid(alpha=0.25)
    figure.suptitle("256×96×288×2 Earth magnetosphere: MHD and x=0, y=20 diagnostics")
    figure.savefig(products / "highres_6h_mhd_x0_y20_diagnostics.png", dpi=150)
    plt.close(figure)

    figure, axes = plt.subplots(3, 1, figsize=(12.5, 9.5), sharex=True, constrained_layout=True)
    axes[0].plot(time_hours, values("cpcp_north_kv"), label="North")
    axes[0].plot(time_hours, values("cpcp_south_kv"), label="South")
    axes[0].set_ylabel("CPCP [kV]")
    axes[0].legend()
    for hemisphere, color in (("north", "tab:blue"), ("south", "tab:orange")):
        axes[1].plot(time_hours, values(f"fac_{hemisphere}_max_uam2"), color=color,
                     label=f"{hemisphere.capitalize()} max")
        axes[1].plot(time_hours, values(f"fac_{hemisphere}_min_uam2"), color=color, ls="--",
                     label=f"{hemisphere.capitalize()} min")
        axes[2].plot(time_hours, values(f"fac_{hemisphere}_rms_uam2"), color=color,
                     label=hemisphere.capitalize())
    axes[1].set_ylabel(r"FAC extrema [$\mu$A/m$^2$]")
    axes[2].set_ylabel(r"FAC RMS [$\mu$A/m$^2$]")
    axes[2].set_xlabel("Physical time [h]")
    axes[1].legend(ncol=2)
    axes[2].legend()
    for axis in axes:
        axis.axvline(1.0, color="0.4", ls="--", lw=0.8)
        axis.axvline(3.0, color="0.4", ls="--", lw=0.8)
        axis.grid(alpha=0.25)
    figure.suptitle("256×96×288×2 Earth magnetosphere: M-I diagnostics (MLT noon convention fixed)")
    figure.savefig(products / "highres_6h_mi_diagnostics.png", dpi=150)
    plt.close(figure)


def validate_final(record: dict[str, object], mi_path: str) -> dict[str, object]:
    checked = []
    for path in (record["p0"], record["p1"]):
        with h5py.File(path, "r") as handle:
            state = np.asarray(handle["state"])
            if not np.isfinite(state).all():
                raise RuntimeError(f"non-finite final compact state: {path}")
            checked.append({"path": str(path), "shape": list(state.shape), "dtype": str(state.dtype)})
    with h5py.File(mi_path, "r") as handle:
        for hemisphere in ("north", "south"):
            for name in ("fac_parallel_A_m2", "potential_V"):
                if not np.isfinite(handle[hemisphere][name][...]).all():
                    raise RuntimeError(f"non-finite final M-I {hemisphere}/{name}")
    return {"compact_mhd": checked, "mi": mi_path, "all_finite": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, nargs=2, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radial-edges", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    output = args.output
    frames = output / "frames"
    products = output / "products"
    if frames.exists() or products.exists():
        raise RuntimeError(f"refusing to overwrite existing postprocess output: {output}")
    frames.mkdir(parents=True)
    products.mkdir(parents=True)

    mhd = discover_mhd(args.inputs)
    mi = discover_mi(args.inputs)
    radial_document = json.loads(args.radial_edges.read_text())
    radial_edges = np.asarray(radial_document["radial_edges"], dtype=float)
    if radial_edges.shape != (257,) or not np.all(np.diff(radial_edges) > 0.0):
        raise RuntimeError(f"invalid radial edges: {radial_edges.shape}")

    tasks = []
    for frame_index, record in enumerate(mhd):
        time_key = int(round(float(record["time_seconds"])))
        tasks.append({
            **record,
            "frame_index": frame_index,
            "mi": mi[time_key],
            "frame": str(frames / f"frame_{frame_index:05d}.png"),
        })

    results = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=initialize_worker,
        initargs=(radial_edges.tolist(),),
    ) as pool:
        for count, result in enumerate(pool.map(render_frame, tasks, chunksize=1), start=1):
            results.append(result)
            if count % 10 == 0 or count == len(tasks):
                print(f"FRAME_PROGRESS {count}/{len(tasks)}", flush=True)
    results.sort(key=lambda record: int(record["frame"]))

    field_order = list(results[0].keys())
    with (products / "highres_6h_diagnostics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=field_order)
        writer.writeheader()
        writer.writerows(results)
    plot_diagnostics(results, products)
    shutil.copy2(results[-1]["frame_path"], products / "highres_6h_final_snapshot.png")

    validation = validate_final(mhd[-1], mi[int(EXPECTED_TIMES[-1])])
    manifest = {
        "status": "frames_complete",
        "grid_per_patch": [256, 96, 288],
        "patches": 2,
        "mpi_ranks": 2048,
        "timeline": {
            "frames": len(results),
            "start_seconds": float(results[0]["time_seconds"]),
            "end_seconds": float(results[-1]["time_seconds"]),
            "cadence_seconds": 120.0,
            "restart_deduplicated": True,
        },
        "inputs": [str(path) for path in args.inputs],
        "fixed_scales": {
            "residual_bz_nt": [-20.0, 20.0],
            "pressure_npa": [0.01, 10.0],
            "fac_uam2": [-4.5, 4.5],
            "potential_contours_kv": [-70.0, 70.0],
        },
        "ionosphere_orientation": "MLT 12 noon at top",
        "validation": validation,
        "final_metrics": results[-1],
    }
    (products / "highres_6h_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
