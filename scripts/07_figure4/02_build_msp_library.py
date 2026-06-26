#!/usr/bin/env python3
"""Build one Paper 2 MSP library from a 0.75-threshold RT input CSV."""

import argparse
import csv
import io
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


BASES = (
    "https://metabolomics-usi.gnps2.org",
    "https://metabolomics-usi.ucsd.edu",
)
ENDPOINTS = (
    ("json", "/json/?usi1={usi}"),
    ("csv", "/csv/?usi1={usi}"),
    ("json", "/json/?usi={usi}"),
    ("json", "/spectrum/json/?usi={usi}"),
)


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value.lower() in {"nan", "none", "na"} else value


def number(value):
    try:
        result = float(value)
        return None if math.isnan(result) else result
    except (TypeError, ValueError):
        return None


def parse_peaks(payload):
    if isinstance(payload, dict):
        peaks = payload.get("peaks")
        if not isinstance(peaks, list):
            spectrum = payload.get("spectrum")
            peaks = spectrum.get("peaks") if isinstance(spectrum, dict) else None
        if isinstance(peaks, list):
            parsed = []
            for peak in peaks:
                if isinstance(peak, (list, tuple)) and len(peak) >= 2:
                    mz, intensity = number(peak[0]), number(peak[1])
                    if mz is not None and intensity is not None:
                        parsed.append((mz, intensity))
            if parsed:
                return parsed
    return None


def parse_peaks_csv(text):
    rows = list(csv.reader(io.StringIO(text)))
    if not rows:
        return None
    start = 0
    if len(rows[0]) >= 2 and (number(rows[0][0]) is None or number(rows[0][1]) is None):
        start = 1
    parsed = []
    for row in rows[start:]:
        if len(row) >= 2:
            mz, intensity = number(row[0]), number(row[1])
            if mz is not None and intensity is not None:
                parsed.append((mz, intensity))
    return parsed or None


def fetch_peaks(usi, timeout, retries):
    encoded = quote(usi, safe="")
    last_error = "no usable peaks"
    for attempt in range(retries):
        for base in BASES:
            for response_type, endpoint in ENDPOINTS:
                url = base + endpoint.format(usi=encoded)
                try:
                    request = Request(url, headers={"User-Agent": "paper2-msp-builder/1.0"})
                    with urlopen(request, timeout=timeout) as response:
                        text = response.read().decode("utf-8")
                    peaks = (
                        parse_peaks(json.loads(text))
                        if response_type == "json"
                        else parse_peaks_csv(text)
                    )
                    if peaks:
                        return peaks, url
                    last_error = f"no usable peaks from {url}"
                except Exception as exc:
                    last_error = f"{type(exc).__name__}: {exc}"
        if attempt + 1 < retries:
            time.sleep(0.5 * (2**attempt))
    raise RuntimeError(last_error)


def normalized_peaks(peaks):
    maximum = max(intensity for _, intensity in peaks) or 1.0
    return sorted(
        ((mz, intensity * 100.0 / maximum) for mz, intensity in peaks),
        key=lambda pair: pair[0],
    )


def msp_block(row, peaks, ion_mode):
    charge = 1 if ion_mode == "P" else -1
    fields = [
        ("Name", clean(row.get("compound.name"))),
        ("PrecursorMZ", number(row.get("mz"))),
        ("RT", number(row.get("rt_for_lib_t075"))),
        ("ExactMass", number(row.get("Monoisotopic.Mass"))),
        ("Spectrum_type", 2),
        ("charge", charge),
        ("Formula", clean(row.get("Formula"))),
        ("SMILES", clean(row.get("smiles"))),
        ("Method", "HCD"),
        ("Collision_energy", "[28.333334]"),
        ("Isolation_window", "4.0"),
        ("Ion_mode", ion_mode),
    ]
    lines = [f"{name}: {value}" for name, value in fields if value != "" and value is not None]
    peaks = normalized_peaks(peaks)
    lines.append(f"Num Peaks: {len(peaks)}")
    lines.extend(f"{mz:.4f}\t{intensity:.1f}" for mz, intensity in peaks)
    return "\n".join(lines) + "\n\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ion-mode", required=True, choices=("P", "N"))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--retry-failures-only", action="store_true")
    args = parser.parse_args()

    started = time.perf_counter()
    with args.input.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"feature.usi", "compound.name", "mz", "rt_for_lib_t075", "Monoisotopic.Mass"}
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    unique_rows = {}
    for row in rows:
        usi = clean(row.get("feature.usi"))
        if usi and usi not in unique_rows:
            unique_rows[usi] = row

    failure_path = args.output.with_name(args.output.stem + "-failures.csv")
    previous_report = None
    if args.retry_failures_only:
        if not failure_path.exists() or not args.output.exists() or not args.report.exists():
            raise SystemExit("Retry mode requires an existing MSP, failure CSV, and report JSON.")
        with failure_path.open(newline="", encoding="utf-8-sig") as handle:
            failed_usis = {clean(row.get("feature.usi")) for row in csv.DictReader(handle)}
        unique_rows = {usi: row for usi, row in unique_rows.items() if usi in failed_usis}
        previous_report = json.loads(args.report.read_text(encoding="utf-8"))

    results, failures = {}, {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(fetch_peaks, usi, args.timeout, args.retries): usi
            for usi in unique_rows
        }
        completed = 0
        for future in as_completed(futures):
            usi = futures[future]
            try:
                results[usi] = future.result()
            except Exception as exc:
                failures[usi] = str(exc)
            completed += 1
            if completed % 250 == 0 or completed == len(futures):
                print(f"[{args.input.name}] fetched {completed}/{len(futures)}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_mode = "a" if args.retry_failures_only else "w"
    with args.output.open(output_mode, encoding="utf-8", newline="\n") as handle:
        for usi, row in unique_rows.items():
            if usi in results:
                peaks, _ = results[usi]
                handle.write(msp_block(row, peaks, args.ion_mode))

    if failures:
        with failure_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=("feature.usi", "error"))
            writer.writeheader()
            writer.writerows(
                {"feature.usi": usi, "error": error} for usi, error in failures.items()
            )
    elif failure_path.exists():
        failure_path.unlink()

    elapsed = time.perf_counter() - started
    previous_entries = previous_report["entries_written"] if previous_report else 0
    previous_elapsed = previous_report["elapsed_seconds"] if previous_report else 0
    report = {
        "input": str(args.input),
        "output": str(args.output),
        "ion_mode": args.ion_mode,
        "input_rows": previous_report["input_rows"] if previous_report else len(rows),
        "unique_usable_usi": previous_report["unique_usable_usi"] if previous_report else len(unique_rows),
        "entries_written": previous_entries + len(results),
        "failures": len(failures),
        "elapsed_seconds": round(previous_elapsed + elapsed, 3),
        "elapsed_minutes": round((previous_elapsed + elapsed) / 60, 3),
        "last_attempt_seconds": round(elapsed, 3),
        "retry_failures_only": args.retry_failures_only,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
