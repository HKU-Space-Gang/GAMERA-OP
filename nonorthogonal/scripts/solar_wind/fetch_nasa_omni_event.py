#!/usr/bin/env python3
"""Download NASA OMNI 1-minute data and create a GAMERA-OP wind HDF5.

The event workflow is intentionally opinionated.  It uses the definitive
OMNI_HRO2_1MIN product, whose timestamps have already been shifted to the
terrestrial bow-shock nose, and writes time relative to simulation start.
Consequently the generated YAML fragment uses a zero time offset and a
reference position of (0, 0, 0).  Users do not calculate or apply a second
L1-to-Earth shift.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Dict, Mapping, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import h5py
import numpy as np


SCRIPT_VERSION = "1.0.0"
NASA_HAPI_DATA = "https://cdaweb.gsfc.nasa.gov/hapi/data"
NASA_DATASET = "OMNI_HRO2_1MIN"
NASA_DOI = "10.48322/mj0k-fq60"
EARTH_RADIUS_KM = 6371.0

# HAPI always returns Time as the first CSV field.  The requested values then
# follow in exactly this order.
NASA_PARAMETERS: Tuple[str, ...] = (
    "BX_GSE",
    "BY_GSM",
    "BZ_GSM",
    "Vx",
    "Vy",
    "Vz",
    "proton_density",
    "T",
)
CSV_COLUMNS: Tuple[str, ...] = ("Time",) + NASA_PARAMETERS

# Values published by the NASA HAPI metadata for OMNI_HRO2_1MIN.  A response
# containing one of these values is missing data, not a physical observation.
NASA_FILL_VALUES: Mapping[str, float] = {
    "BX_GSE": 9999.99,
    "BY_GSM": 9999.99,
    "BZ_GSM": 9999.99,
    "Vx": 99999.9,
    "Vy": 99999.9,
    "Vz": 99999.9,
    "proton_density": 999.99,
    "T": 9999999.0,
}


def parse_utc(text: str) -> datetime:
    normalized = text.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    try:
        value = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ValueError(
            f"invalid UTC timestamp {text!r}; use YYYY-MM-DDTHH:MM:SSZ"
        ) from exc
    if value.tzinfo is None:
        raise ValueError(f"UTC timestamp {text!r} has no timezone")
    return value.astimezone(timezone.utc)


def utc_text(value: datetime) -> str:
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def floor_minute(value: datetime) -> datetime:
    return value.replace(second=0, microsecond=0)


def ceil_minute(value: datetime) -> datetime:
    lower = floor_minute(value)
    return lower if value == lower else lower + timedelta(minutes=1)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_hapi_url(start: datetime, stop_exclusive: datetime) -> str:
    query = urlencode(
        {
            "id": NASA_DATASET,
            "parameters": ",".join(NASA_PARAMETERS),
            "time.min": utc_text(start),
            "time.max": utc_text(stop_exclusive),
            "format": "csv",
        }
    )
    return f"{NASA_HAPI_DATA}?{query}"


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=f".{path.name}.", suffix=".part",
            dir=path.parent, delete=False
        ) as stream:
            temporary_name = stream.name
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary_name, 0o644)
        os.replace(temporary_name, path)
    finally:
        if temporary_name is not None and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def atomic_write_text(path: Path, text: str) -> None:
    atomic_write_bytes(path, text.encode("utf-8"))


def download_once(url: str, path: Path, timeout: float, refresh: bool) -> str:
    """Download one NASA CSV, or reuse the exact local cache."""
    if path.exists() and not refresh:
        return "cached"
    request = Request(
        url,
        headers={
            "Accept": "text/csv",
            "User-Agent": "GAMERA-OP-NASA-OMNI-event-builder/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urlopen(request, timeout=timeout) as response:
                payload = response.read()
            if not payload.strip():
                raise RuntimeError("NASA HAPI returned an empty response")
            if payload.lstrip().startswith(b"{"):
                message = payload.decode("utf-8", errors="replace")[:1000]
                raise RuntimeError(f"NASA HAPI returned JSON instead of CSV: {message}")
            atomic_write_bytes(path, payload)
            return "downloaded"
        except (HTTPError, URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(2.0 ** attempt)
    raise RuntimeError(f"NASA download failed after three attempts: {last_error}")


def is_fill_value(name: str, value: float) -> bool:
    if not math.isfinite(value):
        return True
    fill = NASA_FILL_VALUES[name]
    return math.isclose(value, fill, rel_tol=1.0e-12, abs_tol=1.0e-12)


def parse_hapi_csv(path: Path) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    times = []
    values: Dict[str, list[float]] = {name: [] for name in NASA_PARAMETERS}
    with path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        for line_number, row in enumerate(reader, start=1):
            if not row or (row[0].lstrip().startswith("#")):
                continue
            if row[0].strip().lower() == "time":
                continue
            if len(row) != len(CSV_COLUMNS):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(CSV_COLUMNS)} CSV "
                    f"columns ({','.join(CSV_COLUMNS)}), found {len(row)}"
                )
            sample_time = parse_utc(row[0])
            times.append(sample_time.timestamp())
            for index, name in enumerate(NASA_PARAMETERS, start=1):
                token = row[index].strip()
                try:
                    value = float(token) if token else math.nan
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: {name} is not numeric: {token!r}"
                    ) from exc
                values[name].append(math.nan if is_fill_value(name, value) else value)

    raw_time = np.asarray(times, dtype=np.float64)
    raw_values = {
        name: np.asarray(column, dtype=np.float64)
        for name, column in values.items()
    }
    if raw_time.size == 0:
        raise ValueError(f"NASA CSV contains no data records: {path}")
    if not np.all(np.isfinite(raw_time)) or not np.all(np.diff(raw_time) > 0.0):
        raise ValueError("NASA timestamps are not finite, unique, and increasing")
    for name, column in raw_values.items():
        if column.shape != raw_time.shape:
            raise ValueError(f"NASA column {name} has the wrong length")
    return raw_time, raw_values


def minute_grid(start: datetime, stop: datetime) -> np.ndarray:
    first = int(round(start.timestamp()))
    last = int(round(stop.timestamp()))
    if last < first or (last - first) % 60 != 0:
        raise ValueError("internal minute-grid bounds are invalid")
    return np.arange(first, last + 1, 60, dtype=np.int64).astype(np.float64)


def interpolate_omni(
    raw_time: np.ndarray,
    raw_values: Mapping[str, np.ndarray],
    target_time: np.ndarray,
    max_gap_minutes: int,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float]]:
    cleaned: Dict[str, np.ndarray] = {}
    interpolated: Dict[str, np.ndarray] = {}
    largest_gap: Dict[str, float] = {}

    target_integer = np.rint(target_time).astype(np.int64)
    for name in NASA_PARAMETERS:
        column = raw_values[name]
        valid = np.isfinite(column)
        valid_time = raw_time[valid]
        valid_value = column[valid]
        if valid_time.size < 2:
            raise ValueError(f"NASA variable {name} has fewer than two valid samples")
        if valid_time[0] > target_time[0] or valid_time[-1] < target_time[-1]:
            raise ValueError(
                f"NASA variable {name} does not bracket the requested padded interval"
            )

        separation = np.diff(valid_time)
        overlaps = (valid_time[:-1] <= target_time[-1]) & (
            valid_time[1:] >= target_time[0]
        )
        missing_minutes = np.maximum(separation[overlaps] / 60.0 - 1.0, 0.0)
        maximum = float(np.max(missing_minutes)) if missing_minutes.size else 0.0
        largest_gap[name] = maximum
        if maximum > float(max_gap_minutes) + 1.0e-9:
            raise ValueError(
                f"NASA variable {name} has a {maximum:.0f}-minute missing run; "
                f"the accepted limit is {max_gap_minutes} minutes"
            )

        cleaned[name] = np.interp(target_time, valid_time, valid_value)
        exact_valid_times = set(np.rint(valid_time).astype(np.int64).tolist())
        interpolated[name] = np.asarray(
            [sample not in exact_valid_times for sample in target_integer],
            dtype=np.uint8,
        )
    return cleaned, interpolated, largest_gap


def to_model_native(omni: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    # OMNI velocity is GSE.  OMNI's GSM magnetic Y/Z values are used with the
    # shared GSE/GSM X axis.  The GAMERA-OP native Earth-to-Sun direction is
    # -X, hence only source X components change sign.
    return {
        "D": np.asarray(omni["proton_density"], dtype=np.float64),
        "Temp": np.asarray(omni["T"], dtype=np.float64),
        "Vx": -np.asarray(omni["Vx"], dtype=np.float64),
        "Vy": np.asarray(omni["Vy"], dtype=np.float64),
        "Vz": np.asarray(omni["Vz"], dtype=np.float64),
        "Bx": -np.asarray(omni["BX_GSE"], dtype=np.float64),
        "By": np.asarray(omni["BY_GSM"], dtype=np.float64),
        "Bz": np.asarray(omni["BZ_GSM"], dtype=np.float64),
    }


def validate_model_input(relative_time: np.ndarray,
                         data: Mapping[str, np.ndarray]) -> None:
    if relative_time.size < 2 or not np.all(np.diff(relative_time) > 0.0):
        raise ValueError("generated /T is not a non-empty increasing series")
    if not np.allclose(np.diff(relative_time), 60.0, rtol=0.0, atol=1.0e-9):
        raise ValueError("generated event cadence is not exactly 60 seconds")
    for name, column in data.items():
        if column.shape != relative_time.shape or not np.all(np.isfinite(column)):
            raise ValueError(f"generated {name} is non-finite or has wrong shape")
    if not np.all(data["D"] > 0.0):
        raise ValueError("generated proton density is not strictly positive")
    if not np.all(data["Temp"] > 0.0):
        raise ValueError("generated proton temperature is not strictly positive")
    if not np.all(data["Vx"] > 0.0):
        raise ValueError(
            "OMNI Vx did not map to positive model-native Vx; inspect the raw event"
        )


def companion_paths(output: Path) -> Dict[str, Path]:
    return {
        "yaml": output.with_suffix(".wind.yaml"),
        "qa": output.with_suffix(".qa.png"),
        "receipt": output.with_suffix(".receipt.json"),
    }


def ensure_outputs_available(output: Path, companions: Mapping[str, Path],
                             overwrite: bool, make_plot: bool) -> None:
    paths = [output, companions["yaml"], companions["receipt"]]
    if make_plot:
        paths.append(companions["qa"])
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "refusing to overwrite generated products; use --overwrite or a new "
            f"--output: {', '.join(existing)}"
        )


def write_hdf5(
    output: Path,
    relative_time: np.ndarray,
    source_unix_time: np.ndarray,
    data: Mapping[str, np.ndarray],
    interpolated: Mapping[str, np.ndarray],
    attributes: Mapping[str, object],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as h5:
        h5.create_dataset("T", data=relative_time, dtype=np.float64)
        for name in ("D", "Temp", "Vx", "Vy", "Vz", "Bx", "By", "Bz"):
            h5.create_dataset(name, data=data[name], dtype=np.float64)
        quality = h5.create_group("quality")
        quality.create_dataset(
            "source_unix_seconds", data=source_unix_time, dtype=np.float64
        )
        any_interpolated = np.zeros(relative_time.shape, dtype=np.uint8)
        for source_name in NASA_PARAMETERS:
            mask = interpolated[source_name]
            quality.create_dataset(
                f"interpolated_{source_name}", data=mask, dtype=np.uint8
            )
            any_interpolated = np.maximum(any_interpolated, mask)
        quality.create_dataset(
            "interpolated_any", data=any_interpolated, dtype=np.uint8
        )
        for name, value in attributes.items():
            h5.attrs[name] = value


def write_yaml(path: Path, hdf5_name: str) -> None:
    text = f"""# Generated by fetch_nasa_omni_event.py {SCRIPT_VERSION}
# Copy these keys into the simulation config kept beside {hdf5_name}.
wind_file: {hdf5_name}
wind_input_units: physical
wind_velocity_units: km/s
wind_interpolation: linear
# OMNI_HRO2_1MIN is already shifted to the bow-shock nose.
wind_reference_x: 0.0
wind_reference_y: 0.0
wind_reference_z: 0.0
# /T is written relative to simulation start, so no manual shift is needed.
wind_time_offset: 0.0
wind_enforce_bx_relation: 0
"""
    atomic_write_text(path, text)


def write_qa_plot(
    path: Path,
    relative_time: np.ndarray,
    data: Mapping[str, np.ndarray],
    simulation_duration_seconds: float,
    interpolated_any: np.ndarray,
) -> None:
    try:
        import matplotlib
        # This command-line tool must also work on headless cluster/login
        # nodes.  Select the non-interactive renderer before pyplot import.
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Matplotlib is required for the automatic QA plot; install it or "
            "use --no-qa-plot"
        ) from exc

    hours = relative_time / 3600.0
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True,
                             constrained_layout=True)
    axes[0].plot(hours, data["Bx"], label="native Bx")
    axes[0].plot(hours, data["By"], label="native By")
    axes[0].plot(hours, data["Bz"], label="native Bz")
    axes[0].set_ylabel("B [nT]")
    axes[0].legend(ncol=3)
    axes[1].plot(hours, data["Vx"], label="native Vx")
    axes[1].plot(hours, data["Vy"], label="native Vy")
    axes[1].plot(hours, data["Vz"], label="native Vz")
    axes[1].set_ylabel("V [km/s]")
    axes[1].legend(ncol=3)
    axes[2].plot(hours, data["D"], color="tab:green")
    axes[2].set_ylabel("proton n [cm$^{-3}$]")
    axes[3].plot(hours, data["Temp"], color="tab:red")
    axes[3].set_ylabel("proton T [K]")
    axes[3].set_xlabel("hours relative to simulation start")

    run_end_hours = simulation_duration_seconds / 3600.0
    for axis in axes:
        axis.axvspan(0.0, run_end_hours, color="tab:blue", alpha=0.06)
        axis.grid(alpha=0.25)
    marked = np.flatnonzero(interpolated_any)
    if marked.size:
        axes[3].scatter(hours[marked], data["Temp"][marked], s=8,
                        color="black", label="one or more interpolated fields")
        axes[3].legend(loc="best")
    fig.suptitle("NASA OMNI 1-minute input after GAMERA-OP conversion")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def default_raw_path(output: Path, query_start: datetime,
                     query_stop_exclusive: datetime) -> Path:
    name = (
        f"NASA_{NASA_DATASET}_{query_start:%Y%m%dT%H%M}_"
        f"{query_stop_exclusive:%Y%m%dT%H%M}_UTC.csv"
    )
    return output.parent / "raw" / name


def run_pipeline(args: argparse.Namespace) -> Dict[str, Path]:
    simulation_start = parse_utc(args.simulation_start)
    if args.simulation_stop is not None:
        simulation_stop = parse_utc(args.simulation_stop)
    else:
        simulation_stop = simulation_start + timedelta(hours=args.duration_hours)
    if simulation_stop <= simulation_start:
        raise ValueError("simulation stop must be later than simulation start")
    if args.padding_minutes < 0 or args.max_gap_minutes < 0:
        raise ValueError("padding and maximum gap must be nonnegative")

    coverage_start = floor_minute(
        simulation_start - timedelta(minutes=args.padding_minutes)
    )
    coverage_stop = ceil_minute(
        simulation_stop + timedelta(minutes=args.padding_minutes)
    )
    interpolation_margin = timedelta(minutes=args.max_gap_minutes + 1)
    query_start = coverage_start - interpolation_margin
    # HAPI time.max is exclusive.  Add one cadence after the interpolation
    # margin so the final requested target can be bracketed.
    query_stop_exclusive = coverage_stop + interpolation_margin + timedelta(minutes=1)
    query_url = build_hapi_url(query_start, query_stop_exclusive)

    output = Path(args.output).expanduser().resolve()
    companions = companion_paths(output)
    ensure_outputs_available(
        output, companions, args.overwrite, not args.no_qa_plot
    )
    if args.raw_csv is not None:
        raw_path = Path(args.raw_csv).expanduser().resolve()
        if not raw_path.is_file():
            raise FileNotFoundError(f"--raw-csv does not exist: {raw_path}")
        download_status = "user-supplied-cache"
        retrieval_utc = "not-recorded-for-user-supplied-cache"
    else:
        raw_path = default_raw_path(output, query_start, query_stop_exclusive)
        download_status = download_once(
            query_url, raw_path, args.timeout_seconds, args.refresh
        )
        retrieval_utc = utc_text(datetime.now(timezone.utc))
        source_receipt = raw_path.with_suffix(raw_path.suffix + ".source.json")
        if download_status == "downloaded" or not source_receipt.exists():
            atomic_write_text(
                source_receipt,
                json.dumps(
                    {
                        "dataset": NASA_DATASET,
                        "doi": NASA_DOI,
                        "parameters": list(NASA_PARAMETERS),
                        "query_url": query_url,
                        "retrieved_utc": retrieval_utc,
                        "sha256": sha256_file(raw_path),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
            )
        else:
            try:
                retrieval_utc = json.loads(
                    source_receipt.read_text(encoding="utf-8")
                ).get("retrieved_utc", retrieval_utc)
            except (OSError, json.JSONDecodeError):
                retrieval_utc = "unknown-cached-retrieval-time"

    raw_time, raw_values = parse_hapi_csv(raw_path)
    target_time = minute_grid(coverage_start, coverage_stop)
    omni, interpolated, largest_gap = interpolate_omni(
        raw_time, raw_values, target_time, args.max_gap_minutes
    )
    model_data = to_model_native(omni)
    relative_time = target_time - simulation_start.timestamp()
    validate_model_input(relative_time, model_data)
    minimum_native_vx = float(np.min(model_data["Vx"]))
    maximum_domain_delay_minutes = (
        args.domain_radius_re * EARTH_RADIUS_KM / minimum_native_vx / 60.0
    )
    if args.padding_minutes + 1.0e-9 < maximum_domain_delay_minutes:
        raise ValueError(
            f"{args.padding_minutes} minutes of padding is shorter than the "
            f"{maximum_domain_delay_minutes:.1f}-minute worst-case delay "
            f"audited from native Vx over {args.domain_radius_re:g} RE; "
            "increase --padding-minutes"
        )

    attributes: Dict[str, object] = {
        "builder": "GAMERA-OP fetch_nasa_omni_event.py",
        "builder_version": SCRIPT_VERSION,
        "source_agency": "NASA GSFC SPDF",
        "source_product": NASA_DATASET,
        "source_doi": NASA_DOI,
        "source_hapi_url": query_url,
        "source_parameters": ",".join(NASA_PARAMETERS),
        "source_raw_csv": str(raw_path),
        "source_raw_sha256": sha256_file(raw_path),
        "source_retrieved_utc": retrieval_utc,
        "source_time_status": "OMNI timestamps already shifted to bow-shock nose",
        "simulation_start_utc": utc_text(simulation_start),
        "simulation_stop_utc": utc_text(simulation_stop),
        "event_zero_utc": utc_text(simulation_start),
        "coverage_start_utc": utc_text(coverage_start),
        "coverage_stop_utc": utc_text(coverage_stop),
        "padding_minutes": float(args.padding_minutes),
        "domain_radius_re_for_padding_audit": float(args.domain_radius_re),
        "minimum_model_native_vx_km_s": minimum_native_vx,
        "maximum_domain_delay_minutes": maximum_domain_delay_minutes,
        "accepted_maximum_missing_run_minutes": int(args.max_gap_minutes),
        "gap_policy": "linear interpolation up to accepted limit; longer gaps rejected",
        "coordinate_mapping": "native X=-source X; source Y/Z retained",
        "magnetic_source_frame": "Bx GSE; By/Bz GSM (GSE and GSM share X axis)",
        "velocity_source_frame": "GSE",
        "velocity_unit": "km/s",
        "magnetic_unit": "nT",
        "density_unit": "proton cm^-3",
        "temperature_unit": "K",
        "flow_pressure_used": 0,
        "wind_reference_native_RE": np.asarray([0.0, 0.0, 0.0]),
        "wind_time_offset_code": 0.0,
    }
    write_hdf5(
        output, relative_time, target_time, model_data, interpolated, attributes
    )
    write_yaml(companions["yaml"], output.name)

    interpolated_any = np.zeros(relative_time.shape, dtype=np.uint8)
    for mask in interpolated.values():
        interpolated_any = np.maximum(interpolated_any, mask)
    if not args.no_qa_plot:
        write_qa_plot(
            companions["qa"], relative_time, model_data,
            (simulation_stop - simulation_start).total_seconds(),
            interpolated_any,
        )

    receipt = {
        "builder_version": SCRIPT_VERSION,
        "dataset": NASA_DATASET,
        "dataset_doi": NASA_DOI,
        "download_status": download_status,
        "query_url": query_url,
        "simulation_start_utc": utc_text(simulation_start),
        "simulation_stop_utc": utc_text(simulation_stop),
        "coverage_start_utc": utc_text(coverage_start),
        "coverage_stop_utc": utc_text(coverage_stop),
        "raw_csv": str(raw_path),
        "raw_csv_sha256": sha256_file(raw_path),
        "hdf5": str(output),
        "hdf5_sha256": sha256_file(output),
        "yaml_fragment": str(companions["yaml"]),
        "qa_plot": None if args.no_qa_plot else str(companions["qa"]),
        "sample_count": int(relative_time.size),
        "interpolated_sample_count": int(np.count_nonzero(interpolated_any)),
        "largest_missing_run_minutes_by_variable": largest_gap,
        "padding_audit": {
            "configured_padding_minutes_each_side": int(args.padding_minutes),
            "domain_radius_re": float(args.domain_radius_re),
            "minimum_model_native_vx_km_s": minimum_native_vx,
            "maximum_domain_delay_minutes": maximum_domain_delay_minutes,
            "passed": True,
        },
        "automatic_timing": {
            "event_zero_is_simulation_start": True,
            "omni_already_shifted_to_bow_shock": True,
            "wind_reference_native_RE": [0.0, 0.0, 0.0],
            "wind_time_offset_code": 0.0,
        },
    }
    atomic_write_text(
        companions["receipt"],
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
    )

    print("NASA_OMNI_EVENT_READY")
    print(f"  HDF5:   {output}")
    print(f"  YAML:   {companions['yaml']}")
    print(f"  RAW:    {raw_path} ({download_status})")
    if not args.no_qa_plot:
        print(f"  QA:     {companions['qa']}")
    print(f"  RECEIPT:{companions['receipt']}")
    print(f"  SHA256: {receipt['hdf5_sha256']}")
    return {"hdf5": output, "raw": raw_path, **companions}


def write_self_test_csv(path: Path, start: datetime, stop: datetime) -> None:
    rows = []
    current = start
    index = 0
    while current <= stop:
        temperature = NASA_FILL_VALUES["T"] if index == 7 else 100000.0 + index
        rows.append(
            [
                current.strftime("%Y-%m-%dT%H:%M:%S.000Z"),
                1.0,
                -2.0,
                -5.0,
                -400.0,
                1.0,
                2.0,
                5.0,
                temperature,
            ]
        )
        current += timedelta(minutes=1)
        index += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        csv.writer(stream).writerows(rows)


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="gamera_omni_self_test_") as directory:
        root = Path(directory)
        simulation_start = parse_utc("2015-03-17T00:00:00Z")
        simulation_stop = simulation_start + timedelta(minutes=5)
        coverage_start = floor_minute(simulation_start - timedelta(minutes=2))
        coverage_stop = ceil_minute(simulation_stop + timedelta(minutes=2))
        margin = timedelta(minutes=3)
        raw = root / "fixture.csv"
        write_self_test_csv(raw, coverage_start - margin,
                            coverage_stop + margin + timedelta(minutes=1))
        args = argparse.Namespace(
            simulation_start=utc_text(simulation_start),
            simulation_stop=utc_text(simulation_stop),
            duration_hours=None,
            output=str(root / "event.h5"),
            padding_minutes=2,
            domain_radius_re=1.0,
            max_gap_minutes=2,
            raw_csv=str(raw),
            timeout_seconds=5.0,
            refresh=False,
            overwrite=False,
            no_qa_plot=True,
        )
        products = run_pipeline(args)
        with h5py.File(products["hdf5"], "r") as h5:
            if not np.allclose(h5["Vx"][:], 400.0):
                raise AssertionError("self-test Vx mapping failed")
            if not np.allclose(h5["Bx"][:], -1.0):
                raise AssertionError("self-test Bx mapping failed")
            if int(np.count_nonzero(h5["quality/interpolated_T"][:])) != 1:
                raise AssertionError("self-test gap interpolation audit failed")
            if h5["T"][0] != -120.0 or h5["T"][-1] != 420.0:
                raise AssertionError("self-test automatic timing failed")
        short_padding = argparse.Namespace(**vars(args))
        short_padding.output = str(root / "must_fail_padding.h5")
        short_padding.domain_radius_re = 200.0
        try:
            run_pipeline(short_padding)
        except ValueError as exc:
            if "worst-case delay" not in str(exc):
                raise AssertionError("self-test padding audit failed unexpectedly") from exc
        else:
            raise AssertionError("self-test accepted insufficient event padding")
    print("NASA_OMNI_EVENT_SELF_TEST_PASS")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Download NASA OMNI_HRO2_1MIN once and build a validated "
            "GAMERA-OP physical-units solar-wind HDF5, YAML fragment, QA "
            "plot, raw-data cache, and provenance receipt."
        )
    )
    result.add_argument(
        "--simulation-start",
        help="UTC simulation start, for example 2015-03-17T00:00:00Z",
    )
    stop = result.add_mutually_exclusive_group()
    stop.add_argument("--simulation-stop", help="UTC simulation stop")
    stop.add_argument("--duration-hours", type=float, help="simulation duration")
    result.add_argument("--output", help="output .h5 path")
    result.add_argument(
        "--padding-minutes", type=int, default=120,
        help="automatic OMNI coverage before and after the run (default: 120)",
    )
    result.add_argument(
        "--domain-radius-re", type=float, default=200.0,
        help=(
            "outer domain radius used to audit that padding exceeds the "
            "observed advection delay (default: 200 RE)"
        ),
    )
    result.add_argument(
        "--max-gap-minutes", type=int, default=10,
        help="interpolate missing runs up to this length; reject longer (default: 10)",
    )
    result.add_argument(
        "--raw-csv",
        help="use an existing headerless NASA HAPI CSV instead of downloading",
    )
    result.add_argument(
        "--timeout-seconds", type=float, default=120.0,
        help="NASA request timeout for each of three attempts (default: 120)",
    )
    result.add_argument(
        "--refresh", action="store_true",
        help="download the deterministic raw NASA CSV again instead of using cache",
    )
    result.add_argument(
        "--overwrite", action="store_true",
        help="replace generated HDF5/YAML/QA/receipt products",
    )
    result.add_argument(
        "--no-qa-plot", action="store_true",
        help="skip the automatic Matplotlib QA PNG",
    )
    result.add_argument(
        "--self-test", action="store_true",
        help="run a network-free converter test and exit",
    )
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.self_test:
            self_test()
            return 0
        if not args.simulation_start or not args.output:
            raise ValueError("--simulation-start and --output are required")
        if (args.simulation_stop is None) == (args.duration_hours is None):
            raise ValueError(
                "provide exactly one of --simulation-stop or --duration-hours"
            )
        if args.duration_hours is not None and args.duration_hours <= 0.0:
            raise ValueError("--duration-hours must be positive")
        if args.domain_radius_re <= 0.0:
            raise ValueError("--domain-radius-re must be positive")
        if args.timeout_seconds <= 0.0:
            raise ValueError("--timeout-seconds must be positive")
        run_pipeline(args)
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
