"""Reproduce and verify the 87 main estimates from the released analysis panel.

Run from any working directory with Python 3.14.6 (tested):
    python code/reproduce_main.py

All input paths are relative to this repository. No network access is required.
The calculation matches SI section S7. Coefficients are per one-day earlier
onset; plotting multiplies them by ten. Fraction outcomes remain proportions.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy.stats import t

REPOSITORY = Path(__file__).resolve().parents[1]


def read_panel(manifest_path: Path) -> pd.DataFrame:
    """Verify and concatenate CSV parts without changing their row order."""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    frames = []
    reconstructed = hashlib.sha256()
    header = None
    for index, part in enumerate(manifest["parts"]):
        path = manifest_path.parent / part["file"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != part["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {path.name}")
        with (gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")) as stream:
            first_line = stream.readline()
            if header is None:
                header = first_line
                reconstructed.update(header)
            elif first_line != header:
                raise ValueError(f"Column header differs in {path.name}")
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                reconstructed.update(block)
        frame = pd.read_csv(path, float_precision="round_trip")
        if len(frame) != part["rows"]:
            raise ValueError(f"Row count differs in {path.name}")
        frames.append(frame)
    if reconstructed.hexdigest() != manifest["uncompressed_sha256"]:
        raise ValueError("The concatenated CSV differs from the archived panel")
    panel = pd.concat(frames, ignore_index=True)
    if panel.shape != (manifest["rows"], manifest["columns"]):
        raise ValueError("Unexpected panel dimensions")
    if panel.duplicated(["pixel", "year"]).any():
        raise ValueError("Repeated cell-year identifiers")
    return panel


def group_codes(values: np.ndarray) -> np.ndarray:
    return np.unique(values, return_inverse=True)[1]


def remove_fixed_effects(
    values: np.ndarray,
    groups: tuple[np.ndarray, ...],
    weights: np.ndarray,
    tolerance: float = 1e-8,
    max_iterations: int = 200,
) -> tuple[np.ndarray, int]:
    """Alternately subtract weighted cell and climate-zone–year means."""
    transformed = values.copy()
    for iteration in range(max_iterations):
        previous = transformed.copy()
        for group in groups:
            denominator = np.bincount(group, weights=weights)
            means = np.column_stack(
                [
                    np.bincount(group, weights=weights * transformed[:, j]) / denominator
                    for j in range(transformed.shape[1])
                ]
            )
            transformed -= means[group]
        if np.max(np.abs(transformed - previous)) < tolerance:
            return transformed, iteration + 1
    raise RuntimeError("Fixed-effect projections did not converge")


def cluster_variance(
    residuals: np.ndarray,
    group: np.ndarray,
    weights: np.ndarray,
    predictor: np.ndarray,
    denominator: float,
) -> tuple[float, int]:
    """Cluster-score variance with the archived G/(G−1) correction."""
    codes = group_codes(group)
    count = int(codes.max()) + 1
    if count < 2:
        raise ValueError("At least two clusters are required")
    scores = np.bincount(codes, weights=weights * predictor * residuals)
    variance = np.sum(scores**2) / denominator**2 * count / (count - 1)
    return float(variance), count


def fit_main(panel: pd.DataFrame, outcomes: list[str]) -> tuple[pd.DataFrame, dict]:
    """Fit each outcome on the identical archived common-case sample."""
    required = outcomes + ["pixel", "year", "zone", "lon", "lat", "area_km2", "A"]
    if not np.isfinite(panel[required].to_numpy()).all():
        raise ValueError("The main common-case panel contains missing values")
    if (panel.area_km2 <= 0).any():
        raise ValueError("Area weights must be positive")
    weights = panel.area_km2.to_numpy() / panel.area_km2.mean()
    years = panel.year.to_numpy()
    cells = group_codes(panel.pixel.to_numpy())
    zone_years = group_codes(panel.zone.to_numpy() * 10000 + years)
    matrix = panel[outcomes + ["A"]].to_numpy()
    transformed, iterations = remove_fixed_effects(matrix, (cells, zone_years), weights)
    predictor = transformed[:, -1]
    denominator = np.sum(weights * predictor**2)
    if denominator <= 0:
        raise ValueError("No within-group onset variation")
    estimates = (
        np.sum(weights[:, None] * predictor[:, None] * transformed[:, :-1], axis=0) / denominator
    )
    blocks = np.floor(panel.lon.to_numpy() / 2).astype(int) * 1000 + np.floor(
        panel.lat.to_numpy() / 2
    ).astype(int)
    rows = []
    for j, outcome in enumerate(outcomes):
        residuals = transformed[:, j] - predictor * estimates[j]
        spatial, n_spatial = cluster_variance(residuals, blocks, weights, predictor, denominator)
        temporal, n_years = cluster_variance(residuals, years, weights, predictor, denominator)
        intersection, _ = cluster_variance(
            residuals, blocks * 10000 + years, weights, predictor, denominator
        )
        standard_error = np.sqrt(max(0, spatial + temporal - intersection))
        critical_value = t.ppf(0.975, min(n_spatial, n_years) - 1)
        rows.append(
            {
                "outcome": outcome,
                "beta": estimates[j],
                "ci_low": estimates[j] - critical_value * standard_error,
                "ci_high": estimates[j] + critical_value * standard_error,
                "n": len(panel),
                "pixels": panel.pixel.nunique(),
            }
        )
    details = {
        "iterations": iterations,
        "spatial_clusters": n_spatial,
        "year_clusters": n_years,
        "t_degrees_of_freedom": min(n_spatial, n_years) - 1,
    }
    return pd.DataFrame(rows), details


def verify_results(panel: pd.DataFrame, estimates: pd.DataFrame, reference: pd.DataFrame) -> dict:
    if len(estimates) != 87 or len(reference) != 87:
        raise ValueError("Exactly 87 main outcomes must be compared")
    if set(estimates.outcome) != set(reference.outcome):
        raise ValueError("Main outcome names differ")
    joined = estimates.merge(
        reference, on="outcome", suffixes=("_new", "_reference"), validate="one_to_one"
    )
    error = max(
        np.max(np.abs(joined[f"{key}_new"] - joined[f"{key}_reference"]))
        for key in ("beta", "ci_low", "ci_high")
    )
    if error >= 1e-7:
        raise ValueError(f"Estimate or interval mismatch: {error}")
    for key in ("n", "pixels"):
        if not (joined[f"{key}_new"] == joined[f"{key}_reference"]).all():
            raise ValueError(f"Sample count mismatch: {key}")
    identity_errors = {}
    for suffix in ("I", "K_loo", "R_loo", "early", "late", "annual"):
        difference = panel[f"E_{suffix}"] - panel[f"Et_{suffix}"] - panel[f"En_{suffix}"]
        identity_errors[suffix] = float(np.nanmax(np.abs(difference)))
        if identity_errors[suffix] >= 1e-3:
            raise ValueError(f"E=T+N identity failed: {suffix}")
    return {
        "maximum_estimate_or_interval_difference": float(error),
        "component_identity_errors_mm": identity_errors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, help="Optional complete Dataset S1 CSV or CSV.gz")
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results/reproduced_main.csv")
    args = parser.parse_args()
    panel = (
        pd.read_csv(args.panel, float_precision="round_trip")
        if args.panel
        else read_panel(REPOSITORY / "data/panel/manifest.json")
    )
    reference = pd.read_csv(REPOSITORY / "reference_results/main_coefficients.csv")
    outcome_set = set(reference.outcome)
    # Preserve source-column order, as in the original calculation.
    outcomes = [key for key in panel.columns if key in outcome_set]
    estimates, details = fit_main(panel, outcomes)
    checks = verify_results(panel, estimates, reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    estimates.to_csv(args.output, index=False)
    report = {
        "status": "PASS",
        "outcomes": len(estimates),
        "rows": len(panel),
        "cells": int(panel.pixel.nunique()),
        **details,
        **checks,
    }
    args.output.with_suffix(".verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
