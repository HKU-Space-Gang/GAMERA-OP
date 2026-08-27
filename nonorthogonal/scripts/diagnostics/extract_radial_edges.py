#!/usr/bin/env python3
"""Extract the analytic radial coordinate from the compact patch grid."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("grid", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with h5py.File(args.grid, "r") as handle:
        vertices = np.asarray(handle["vertices"], dtype=float)
        cells = tuple(int(value) for value in handle["global_cells"][...])
    if vertices.shape != (3, cells[0] + 1, cells[1] + 1, cells[2] + 1):
        raise RuntimeError(f"unexpected vertex shape {vertices.shape}")
    radial_edges = np.linalg.norm(vertices[:, :, 0, 0], axis=0)
    if radial_edges.shape != (cells[0] + 1,) or not np.all(np.diff(radial_edges) > 0.0):
        raise RuntimeError("invalid radial coordinate")
    document = {
        "source": str(args.grid),
        "radial_edges": radial_edges.tolist(),
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n")


if __name__ == "__main__":
    main()
