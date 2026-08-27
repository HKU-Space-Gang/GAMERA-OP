#!/usr/bin/env python3
"""Render a full North/South M-I physics snapshot with the native MLT view."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np


def cyclic(field: np.ndarray, longitude: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.concatenate((longitude, [longitude[0] + 2.0 * np.pi])),
        np.concatenate((field, field[:, :1]), axis=1),
    )


def polar_axis(axis: plt.Axes, colatitude: np.ndarray) -> None:
    # Native M-I longitude: phi=0 is 00 MLT and phi=pi is 12 MLT.
    axis.set_theta_zero_location("S")
    axis.set_theta_direction(1)
    axis.set_thetagrids([0, 90, 180, 270], ["00", "06", "12", "18"])
    axis.set_ylim(0.0, float(colatitude[-1]))
    ticks_deg = np.arange(10.0, np.degrees(colatitude[-1]) + 0.1, 10.0)
    axis.set_yticks(np.radians(ticks_deg))
    axis.set_yticklabels([f"{90.0-value:.0f}\N{DEGREE SIGN}" for value in ticks_deg])
    axis.set_rlabel_position(225)
    axis.grid(alpha=0.28, linewidth=0.45)


def draw(
    axis: plt.Axes,
    field: np.ndarray,
    longitude: np.ndarray,
    colatitude: np.ndarray,
    boundary_mlat: np.ndarray,
    title: str,
    cmap: str,
    norm: colors.Normalize,
) -> object:
    lon, values = cyclic(field, longitude)
    ll, cc = np.meshgrid(lon, colatitude)
    image = axis.pcolormesh(ll, cc, values, shading="auto", cmap=cmap, norm=norm)
    boundary_lon = np.concatenate((longitude, [longitude[0] + 2.0 * np.pi]))
    boundary_radius = np.radians(
        90.0 - np.concatenate((boundary_mlat, [boundary_mlat[0]]))
    )
    axis.plot(boundary_lon, boundary_radius, color="#38ff42", lw=1.1)
    # Native terminator for sunward=(-X): dawn/dusk meridians at 06/18 MLT.
    axis.plot([np.pi / 2.0, np.pi / 2.0], [0.0, colatitude[-1]],
              color="0.75", ls="--", lw=0.7)
    axis.plot([3.0 * np.pi / 2.0, 3.0 * np.pi / 2.0],
              [0.0, colatitude[-1]], color="0.75", ls="--", lw=0.7)
    polar_axis(axis, colatitude)
    axis.set_title(title, fontsize=9.5, pad=13)
    return image


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(args.input, "r") as handle:
        time_seconds = float(handle.attrs["time_seconds"])
        longitude = np.asarray(handle["longitude_rad"], dtype=float)
        colatitude = np.asarray(handle["colatitude_rad"], dtype=float)
        sunward = np.asarray(
            [handle.attrs[f"sunward_direction_{axis}"] for axis in "xyz"],
            dtype=float,
        )
        sunward /= np.linalg.norm(sunward)
        data: dict[str, dict[str, np.ndarray]] = {}
        multipliers: dict[str, int] = {}
        for hemisphere, fallback in (("north", -1), ("south", 1)):
            group = handle[hemisphere]
            multipliers[hemisphere] = int(
                group.attrs.get("upward_fac_multiplier", fallback)
            )
            data[hemisphere] = {
                name: np.asarray(group[name], dtype=float)
                for name in (
                    "potential_V", "fac_parallel_A_m2", "Sigma_P_S", "Sigma_H_S",
                    "F_E_erg_cm2_s", "F_N_cm2_s", "DPB_boundary_MLAT_deg",
                )
            }

    potential = {
        hemisphere: values["potential_V"] * 1.0e-3
        for hemisphere, values in data.items()
    }
    fac = {
        hemisphere: values["fac_parallel_A_m2"] * multipliers[hemisphere] * 1.0e6
        for hemisphere, values in data.items()
    }
    p_limit = max(float(np.max(np.abs(values))) for values in potential.values())
    p_limit = max(5.0 * np.ceil(p_limit / 5.0), 5.0)
    maxima = {
        name: max(float(np.max(values[name])) for values in data.values())
        for name in ("Sigma_P_S", "Sigma_H_S", "F_E_erg_cm2_s", "F_N_cm2_s")
    }
    norms = {
        "potential": colors.TwoSlopeNorm(vmin=-p_limit, vcenter=0.0, vmax=p_limit),
        "fac": colors.TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5),
        "Sigma_P_S": colors.Normalize(0.0, max(maxima["Sigma_P_S"], 1.0e-12)),
        "Sigma_H_S": colors.Normalize(0.0, max(maxima["Sigma_H_S"], 1.0e-12)),
        "F_E_erg_cm2_s": colors.Normalize(0.0, max(maxima["F_E_erg_cm2_s"], 1.0e-30)),
        "F_N_cm2_s": colors.Normalize(0.0, max(maxima["F_N_cm2_s"], 1.0e-30)),
    }

    figure, axes = plt.subplots(
        4, 3, figsize=(13.8, 16.4), subplot_kw={"projection": "polar"},
        constrained_layout=True,
    )
    for offset, hemisphere in ((0, "north"), (2, "south")):
        label = hemisphere.capitalize()
        values = data[hemisphere]
        panels = (
            (potential[hemisphere], rf"{label}: potential $\Phi$ [kV]", "RdBu_r", norms["potential"]),
            (fac[hemisphere], rf"{label}: upward FAC [$\mu$A m$^{{-2}}$]", "RdBu_r", norms["fac"]),
            (values["Sigma_P_S"], rf"{label}: $\Sigma_P$ [S]", "viridis", norms["Sigma_P_S"]),
            (values["Sigma_H_S"], rf"{label}: $\Sigma_H$ [S]", "viridis", norms["Sigma_H_S"]),
            (values["F_E_erg_cm2_s"], rf"{label}: $F_E$ [erg cm$^{{-2}}$ s$^{{-1}}$] (linear)", "inferno", norms["F_E_erg_cm2_s"]),
            (values["F_N_cm2_s"], rf"{label}: $F_N$ [cm$^{{-2}}$ s$^{{-1}}$] (linear)", "inferno", norms["F_N_cm2_s"]),
        )
        for local_index, (field, title, cmap, norm) in enumerate(panels):
            row = offset + local_index // 3
            column = local_index % 3
            image = draw(
                axes[row, column], field, longitude, colatitude,
                values["DPB_boundary_MLAT_deg"], title, cmap, norm,
            )
            figure.colorbar(image, ax=axes[row, column], pad=0.075, shrink=0.73)

    cpcp_n = float(np.ptp(potential["north"]))
    cpcp_s = float(np.ptp(potential["south"]))
    figure.suptitle(
        f"M-I physics quick snapshot — t={time_seconds:.0f} s "
        f"({time_seconds / 3600.0:.2f} h), grid={longitude.size}x{colatitude.size}/hemisphere\n"
        f"CPCP N/S={cpcp_n:.1f}/{cpcp_s:.1f} kV; "
        f"native Earth-to-Sun=({sunward[0]:+.1f},{sunward[1]:+.1f},{sunward[2]:+.1f}); "
        "12 noon top, 06 dawn right, 00 midnight bottom, 18 dusk left; green=DPB",
        fontsize=13,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=175, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
