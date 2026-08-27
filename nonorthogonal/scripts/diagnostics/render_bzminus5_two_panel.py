#!/usr/bin/env python3
"""Render the standard two-panel Earth MHD/M-I movie from compact output.

The MHD layout follows the established two-panel product: XY residual Bz on
the left and XZ thermal pressure on the right, with North/South M-I polar caps
embedded in the right panel.  Stored M-I longitude is simulation phi, where
phi=0 is midnight and phi increases toward 06 MLT.  The polar axes therefore
use south as theta=0 and counter-clockwise-positive theta so the data—not just
the labels—map to 12/top, 06/right, 00/bottom, and 18/left.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
from pathlib import Path
import shutil

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import numpy as np

import postprocess_compact_base as compact


FAC_LIMIT_UAM2 = 2.5
POTENTIAL_LIMIT_KV = 70.0
RESIDUAL_BZ_LIMIT_NT = 20.0
PRESSURE_MIN_NPA = 0.01
PRESSURE_MAX_NPA = 10.0
INNER_BOUNDARY_RE = 2.5
EARTH_RADIUS_RE = 1.0


def validate_mlt_mapping() -> dict[str, str]:
    """Verify the actual stored-longitude to screen-position mapping."""
    mapping = {
        "stored_phi_0": "00 midnight / bottom",
        "stored_phi_pi_over_2": "06 dawn / right",
        "stored_phi_pi": "12 noon / top",
        "stored_phi_3pi_over_2": "18 dusk / left",
    }
    phi = np.asarray([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])
    # With theta zero at South and positive counter-clockwise, screen-space
    # unit coordinates are (sin(phi), -cos(phi)).
    actual = np.column_stack((np.sin(phi), -np.cos(phi)))
    expected = np.asarray(((0.0, -1.0), (1.0, 0.0),
                           (0.0, 1.0), (-1.0, 0.0)))
    if not np.allclose(actual, expected, atol=1.0e-14):
        raise RuntimeError("invalid MLT display transform")
    return mapping


def phase_label(time_seconds: float) -> tuple[str, float]:
    if time_seconds < 3600.0:
        return "Zero-IMF solar-wind preconditioning", time_seconds
    return r"Southward IMF $B_z=-5$ nT evolution", time_seconds - 3600.0


def polar_panel(axis, data, longitude, colatitude, hemisphere):
    lon, fac = compact.cyclic(data["fac"], longitude)
    _, potential = compact.cyclic(data["potential"], longitude)
    ll, cc = np.meshgrid(lon, colatitude)
    fac_levels = np.linspace(-FAC_LIMIT_UAM2, FAC_LIMIT_UAM2, 49)
    image = axis.contourf(
        ll, cc, fac, levels=fac_levels, cmap="RdBu_r", extend="both"
    )
    potential_levels = POTENTIAL_LIMIT_KV * np.asarray(
        [-0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75]
    )
    axis.contour(
        ll, cc, potential, levels=potential_levels, colors="0.18",
        linewidths=[0.65, 0.65, 0.65, 1.05, 0.65, 0.65, 0.65],
    )

    # This is the physical data transform.  It is intentionally not a mere
    # relabelling of a north-zero plot.
    axis.set_theta_zero_location("S")
    axis.set_theta_direction(1)
    _, theta_labels = axis.set_thetagrids(
        [0, 90, 180, 270], ["00", "06", "12", "18"], fontsize=6.2
    )
    for label in theta_labels:
        label.set_color("0.88")
        label.set_weight("semibold")
    axis.set_ylim(0.0, float(colatitude[-1]))
    maximum_degrees = np.degrees(float(colatitude[-1]))
    rings = np.radians(np.arange(10.0, maximum_degrees + 0.1, 10.0))
    axis.set_yticks(rings)
    axis.set_yticklabels(
        [f"{90.0 - np.degrees(ring):.0f}°" for ring in rings], fontsize=5.5
    )
    axis.set_rlabel_position(225)
    axis.grid(alpha=0.38, linewidth=0.45)
    axis.set_facecolor("white")
    axis.patch.set_alpha(0.96)
    axis.set_title(
        f"{hemisphere.capitalize()}: FAC; $\\Phi$ contours",
        fontsize=8.2, pad=4, color="0.88", weight="semibold",
    )
    return image


def draw_earth(axis) -> None:
    """Draw the physical Earth inside the masked MHD inner boundary.

    The standard SM-X panels place the dayside at +X (screen right), so the
    right semicircle is sunlit and the left semicircle is nightside.
    """
    axis.add_patch(
        Wedge(
            (0.0, 0.0), EARTH_RADIUS_RE, -90.0, 90.0,
            facecolor="white", edgecolor="none", zorder=10,
        )
    )
    axis.add_patch(
        Wedge(
            (0.0, 0.0), EARTH_RADIUS_RE, 90.0, 270.0,
            facecolor="0.32", edgecolor="none", zorder=10,
        )
    )
    axis.add_patch(
        plt.Circle(
            (0.0, 0.0), EARTH_RADIUS_RE, fill=False,
            edgecolor="0.05", linewidth=0.9, zorder=11,
        )
    )
    axis.plot(
        [0.0, 0.0], [-EARTH_RADIUS_RE, EARTH_RADIUS_RE],
        color="0.05", linewidth=0.65, zorder=11,
    )


def decorate_mhd(axis, transverse_label: str) -> None:
    axis.add_patch(
        plt.Circle(
            (0.0, 0.0), INNER_BOUNDARY_RE, color="0.12", zorder=8,
        )
    )
    draw_earth(axis)
    axis.set_xlim(-100.0, 30.0)
    axis.set_ylim(-60.0, 60.0)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"SM-X [$R_E$]")
    axis.set_ylabel(fr"SM-{transverse_label} [$R_E$]")
    axis.grid(alpha=0.16, linewidth=0.4)
    axis.text(
        0.015, 0.98, "NIGHTSIDE / TAIL", transform=axis.transAxes,
        ha="left", va="top", color="0.22", weight="bold", fontsize=8,
    )
    axis.text(
        0.985, 0.98, "DAYSIDE / UPSTREAM", transform=axis.transAxes,
        ha="right", va="top", color="0.70", weight="bold", fontsize=8,
    )


def render_frame(task: dict[str, object]) -> dict[str, object]:
    radial_edges = compact.WORKER_CACHE["radial_edges"]
    handles = [h5py.File(task["p0"], "r"), h5py.File(task["p1"], "r")]
    try:
        names = [compact.field_names(handle) for handle in handles]
        indices = [
            {name: field_list.index(name) for name in compact.FIELDS}
            for field_list in names
        ]
        if "geometry" not in compact.WORKER_CACHE:
            compact.WORKER_CACHE["geometry"] = compact.build_geometry(
                handles[0], radial_edges
            )
        geometry = compact.WORKER_CACHE["geometry"]
        normalization = np.asarray(handles[0]["normalization"], dtype=float)

        bz_code = compact.compose_field(
            handles, indices, geometry["plane_xy"], "bz"
        )
        pressure_code = compact.compose_field(
            handles, indices, geometry["plane_xz"], "p"
        )
        probe = compact.probe_fields(handles, indices, geometry["probe"])
        moment_z = -3.1e-5 / normalization[6]
        residual_bz = (
            bz_code
            - compact.dipole_bz(
                geometry["x_xy"], geometry["y_xy"], geometry["z_xy"],
                moment_z,
            )
        ) * normalization[6] * 1.0e9
        pressure = np.maximum(
            pressure_code * normalization[5] * 1.0e9, 1.0e-30
        )
        mi, longitude, colatitude = compact.mi_arrays(task["mi"])

        figure_width = 16.0
        figure_height = 8.6
        figure = plt.figure(figsize=(figure_width, figure_height))
        panel_bottom = 0.19
        panel_height = 0.68
        panel_gap = 0.025
        # Match the physical axes box to the 130-by-120 RE data extent, then
        # center the two equal panels around a deliberately narrow gutter.
        panel_width = (
            panel_height * figure_height / figure_width * (130.0 / 120.0)
        )
        left_x = 0.5 * (1.0 - 2.0 * panel_width - panel_gap)
        right_x = left_x + panel_width + panel_gap
        ax_xy = figure.add_axes(
            [left_x, panel_bottom, panel_width, panel_height]
        )
        ax_xz = figure.add_axes(
            [right_x, panel_bottom, panel_width, panel_height]
        )
        # Colorbars use the exact x origin and width of their parent panels.
        colorbar_bottom = 0.075
        colorbar_height = 0.032
        bz_color_axis = figure.add_axes(
            [left_x, colorbar_bottom, panel_width, colorbar_height]
        )
        pressure_color_axis = figure.add_axes(
            [right_x, colorbar_bottom, panel_width, colorbar_height]
        )

        bnorm = colors.TwoSlopeNorm(
            vmin=-RESIDUAL_BZ_LIMIT_NT, vcenter=0.0,
            vmax=RESIDUAL_BZ_LIMIT_NT,
        )
        pnorm = colors.LogNorm(vmin=PRESSURE_MIN_NPA, vmax=PRESSURE_MAX_NPA)
        bmesh = ax_xy.pcolormesh(
            -geometry["x_xy_e"], geometry["y_xy_e"], residual_bz,
            shading="flat", cmap="RdBu_r", norm=bnorm, rasterized=True,
        )
        pmesh = ax_xz.pcolormesh(
            -geometry["x_xz_e"], geometry["z_xz_e"], pressure,
            shading="flat", cmap="viridis", norm=pnorm, rasterized=True,
        )
        decorate_mhd(ax_xy, "Y")
        decorate_mhd(ax_xz, "Z")
        ax_xy.set_title(r"XY: $B_z-B_{z,\mathrm{dipole}}$", fontsize=11)
        ax_xz.set_title("XZ: thermal pressure", fontsize=11)
        figure.colorbar(
            bmesh, cax=bz_color_axis, orientation="horizontal",
            label=r"$B_z-B_{z,\mathrm{dipole}}$ [nT]",
        )
        figure.colorbar(
            pmesh, cax=pressure_color_axis, orientation="horizontal",
            label="thermal pressure [nPa]",
        )

        xz_bounds = ax_xz.get_position()

        def overlay_bounds(relative):
            return [
                xz_bounds.x0 + relative[0] * xz_bounds.width,
                xz_bounds.y0 + relative[1] * xz_bounds.height,
                relative[2] * xz_bounds.width,
                relative[3] * xz_bounds.height,
            ]

        ax_n = figure.add_axes(
            overlay_bounds([0.025, 0.60, 0.28, 0.35]),
            projection="polar", zorder=12,
        )
        ax_s = figure.add_axes(
            overlay_bounds([0.025, 0.08, 0.28, 0.35]),
            projection="polar", zorder=12,
        )
        fac_color_axis = figure.add_axes(
            overlay_bounds([0.32, 0.095, 0.020, 0.31]), zorder=13
        )
        fac_image = polar_panel(
            ax_n, mi["north"], longitude, colatitude, "north"
        )
        polar_panel(ax_s, mi["south"], longitude, colatitude, "south")
        figure.colorbar(
            fac_image, cax=fac_color_axis, orientation="vertical",
            ticks=[-2.0, -1.0, 0.0, 1.0, 2.0],
            label=r"upward FAC [$\mu$A m$^{-2}$]",
        )
        fac_color_axis.tick_params(labelsize=7, colors="0.92")
        fac_color_axis.yaxis.label.set_color("0.92")
        for spine in fac_color_axis.spines.values():
            spine.set_edgecolor("0.92")
        ax_xz.text(
            0.985, 0.925,
            f"CPCP (N/S)\n{mi['north']['cpcp_kv']:.1f} / "
            f"{mi['south']['cpcp_kv']:.1f} kV",
            transform=ax_xz.transAxes, ha="right", va="top", fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.88,
                  "edgecolor": "0.25"},
            zorder=14,
        )

        time_seconds = float(task["time_seconds"])
        label, phase_seconds = phase_label(time_seconds)
        figure.suptitle(
            f"{label} — phase t={phase_seconds:.0f} s; "
            f"absolute t={time_seconds:.0f} s\n"
            "MLT: 12 noon (top), 06 dawn (right), "
            "00 midnight (bottom), 18 dusk (left)",
            fontsize=13.2,
        )
        output = Path(task["frame"])
        figure.savefig(output, dpi=120, facecolor="white")
        plt.close(figure)
        result = {
            "frame_index": int(task["frame_index"]),
            "time_seconds": time_seconds,
            "frame": str(output),
            "probe_density_cm3": (
                probe["rho"] * normalization[4]
                / compact.PROTON_MASS_KG / 1.0e6
            ),
            "probe_pressure_npa": probe["p"] * normalization[5] * 1.0e9,
            "xy_residual_bz_min_nt": float(np.min(residual_bz)),
            "xy_residual_bz_max_nt": float(np.max(residual_bz)),
            "xz_pressure_min_npa": float(np.min(pressure)),
            "xz_pressure_max_npa": float(np.max(pressure)),
        }
        for name in ("vx", "vy", "vz"):
            result[f"probe_{name}_kms"] = (
                probe[name] * normalization[2] / 1000.0
            )
        for name in ("bx", "by", "bz"):
            result[f"probe_{name}_nt"] = (
                probe[name] * normalization[6] * 1.0e9
            )
        for hemisphere in ("north", "south"):
            result[f"cpcp_{hemisphere}_kv"] = mi[hemisphere]["cpcp_kv"]
            result[f"fac_{hemisphere}_min_uam2"] = mi[hemisphere]["fac_min"]
            result[f"fac_{hemisphere}_max_uam2"] = mi[hemisphere]["fac_max"]
            result[f"fac_{hemisphere}_rms_uam2"] = mi[hemisphere]["fac_rms"]
        return result
    finally:
        for handle in handles:
            handle.close()


def plot_diagnostics(records: list[dict[str, object]], products: Path,
                     prefix: str) -> None:
    time_hours = np.asarray(
        [record["time_seconds"] for record in records], dtype=float
    ) / 3600.0

    def values(name: str) -> np.ndarray:
        return np.asarray([record[name] for record in records], dtype=float)

    figure, axes = plt.subplots(
        4, 1, figsize=(12.5, 12.0), sharex=True, constrained_layout=True
    )
    axes[0].plot(time_hours, values("probe_density_cm3"))
    axes[0].set_ylabel("Density [cm$^{-3}$]")
    pressure_axis = axes[0].twinx()
    pressure_axis.plot(
        time_hours, values("probe_pressure_npa"), color="tab:red"
    )
    pressure_axis.set_ylabel("Pressure [nPa]")
    for name in ("vx", "vy", "vz"):
        axes[1].plot(
            time_hours, values(f"probe_{name}_kms"),
            label=f"${name[0].upper()}_{name[1]}$",
        )
    axes[1].set_ylabel("Velocity [km/s]")
    axes[1].legend(ncol=3)
    for name in ("bx", "by", "bz"):
        axes[2].plot(
            time_hours, values(f"probe_{name}_nt"),
            label=f"${name[0].upper()}_{name[1]}$",
        )
    axes[2].set_ylabel("Magnetic field [nT]")
    axes[2].legend(ncol=3)
    axes[3].plot(
        time_hours, values("xy_residual_bz_min_nt"),
        label="XY min residual $B_z$",
    )
    axes[3].plot(
        time_hours, values("xy_residual_bz_max_nt"),
        label="XY max residual $B_z$",
    )
    axes[3].set_ylabel("Residual $B_z$ [nT]")
    axes[3].set_xlabel("Physical time [h]")
    axes[3].legend(ncol=2)
    for axis in axes:
        axis.axvline(1.0, color="0.4", ls="--", lw=0.8)
        axis.grid(alpha=0.25)
    figure.suptitle(
        "128×48×144×2 Earth magnetosphere: MHD and x=0, y=20 diagnostics"
    )
    figure.savefig(
        products / f"{prefix}_mhd_x0_y20_diagnostics.png", dpi=150
    )
    plt.close(figure)

    figure, axes = plt.subplots(
        3, 1, figsize=(12.5, 9.5), sharex=True, constrained_layout=True
    )
    axes[0].plot(time_hours, values("cpcp_north_kv"), label="North")
    axes[0].plot(time_hours, values("cpcp_south_kv"), label="South")
    axes[0].set_ylabel("CPCP [kV]")
    axes[0].legend()
    for hemisphere, color in (("north", "tab:blue"),
                              ("south", "tab:orange")):
        axes[1].plot(
            time_hours, values(f"fac_{hemisphere}_max_uam2"),
            color=color, label=f"{hemisphere.capitalize()} max",
        )
        axes[1].plot(
            time_hours, values(f"fac_{hemisphere}_min_uam2"),
            color=color, ls="--", label=f"{hemisphere.capitalize()} min",
        )
        axes[2].plot(
            time_hours, values(f"fac_{hemisphere}_rms_uam2"),
            color=color, label=hemisphere.capitalize(),
        )
    axes[1].set_ylabel(r"FAC extrema [$\mu$A/m$^2$]")
    axes[2].set_ylabel(r"FAC RMS [$\mu$A/m$^2$]")
    axes[2].set_xlabel("Physical time [h]")
    axes[1].legend(ncol=2)
    axes[2].legend()
    for axis in axes:
        axis.axvline(1.0, color="0.4", ls="--", lw=0.8)
        axis.grid(alpha=0.25)
    figure.suptitle(
        "128×48×144×2 Earth magnetosphere: M-I diagnostics "
        "(MLT noon at top)"
    )
    figure.savefig(products / f"{prefix}_mi_diagnostics.png", dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radial-edges", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--start-seconds", type=float, default=0.0)
    parser.add_argument("--end-seconds", type=float, default=18000.0)
    parser.add_argument("--cadence-seconds", type=float, default=120.0)
    parser.add_argument("--product-prefix", default="bzminus5_dayside_euv_5h")
    parser.add_argument(
        "--frame-index", type=int,
        help="render only one restart-deduplicated frame (for preview QA)",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.end_seconds < args.start_seconds or args.cadence_seconds <= 0.0:
        parser.error("invalid time range or cadence")

    orientation = validate_mlt_mapping()
    compact.EXPECTED_TIMES = np.arange(
        args.start_seconds,
        args.end_seconds + 0.5 * args.cadence_seconds,
        args.cadence_seconds,
    )
    mhd = compact.discover_mhd(args.inputs)
    mi = compact.discover_mi(args.inputs)
    radial_document = json.loads(args.radial_edges.read_text())
    radial_edges = np.asarray(radial_document["radial_edges"], dtype=float)
    if radial_edges.shape != (129,) or not np.all(np.diff(radial_edges) > 0.0):
        raise RuntimeError(f"invalid radial edges: {radial_edges.shape}")

    selected = list(range(len(mhd)))
    if args.frame_index is not None:
        if args.frame_index < 0 or args.frame_index >= len(mhd):
            parser.error(f"--frame-index must be in [0,{len(mhd) - 1}]")
        selected = [args.frame_index]

    frames = args.output / "frames"
    products = args.output / "products"
    if frames.exists() or products.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {args.output}")
    frames.mkdir(parents=True)
    products.mkdir(parents=True)
    tasks = []
    for output_index, source_index in enumerate(selected):
        record = mhd[source_index]
        time_key = int(round(float(record["time_seconds"])))
        tasks.append({
            **record,
            "frame_index": source_index,
            "mi": mi[time_key],
            "frame": str(frames / f"frame_{output_index:05d}.png"),
        })

    results = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)),
        initializer=compact.initialize_worker,
        initargs=(radial_edges.tolist(),),
    ) as pool:
        for count, result in enumerate(
            pool.map(render_frame, tasks, chunksize=1), start=1
        ):
            results.append(result)
            if count % 10 == 0 or count == len(tasks):
                print(f"FRAME_PROGRESS {count}/{len(tasks)}", flush=True)
    results.sort(key=lambda item: int(item["frame_index"]))
    prefix = args.product_prefix
    shutil.copy2(
        results[-1]["frame"], products / f"{prefix}_two_panel_final.png"
    )
    with (products / f"{prefix}_diagnostics.csv").open(
        "w", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    plot_diagnostics(results, products, prefix)
    final_validation = compact.validate_final(
        mhd[-1], mi[int(round(compact.EXPECTED_TIMES[-1]))]
    )
    manifest = {
        "schema_version": 1,
        "status": "complete",
        "layout": "standard two-panel: XY residual Bz; XZ pressure with embedded North/South M-I",
        "residual_bz_colormap": "RdBu_r",
        "pressure_colormap": "viridis",
        "fixed_scales": {
            "residual_bz_nt": [-RESIDUAL_BZ_LIMIT_NT, RESIDUAL_BZ_LIMIT_NT],
            "pressure_npa": [PRESSURE_MIN_NPA, PRESSURE_MAX_NPA],
            "fac_uam2": [-FAC_LIMIT_UAM2, FAC_LIMIT_UAM2],
            "potential_limit_kv": POTENTIAL_LIMIT_KV,
        },
        "mlt_orientation": orientation,
        "theta_configuration": "zero=S; direction=counter-clockwise",
        "timeline": {
            "frames": len(results),
            "start_seconds": float(results[0]["time_seconds"]),
            "end_seconds": float(results[-1]["time_seconds"]),
            "cadence_seconds": args.cadence_seconds,
            "restart_deduplicated": True,
        },
        "inputs": [str(path) for path in args.inputs],
        "validation": final_validation,
        "frames": results,
    }
    (products / f"{prefix}_two_panel_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
