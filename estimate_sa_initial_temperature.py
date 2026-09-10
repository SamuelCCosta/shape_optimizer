#!/usr/bin/env python3
"""Estimate simulated-annealing initial temperatures from tracked runs.

The script reads run metadata from a configurable SQLite table (default: ``results``).
For every run it resolves ``track_file_name``, reconstructs proposed cost
differences from the accepted-state history, and estimates T0 from a chosen
initial window of proposals.

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import random
import sqlite3
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SCRIPT_VERSION = "1.0.0"
DEFAULT_TARGETS = (0.80, 0.85, 0.90, 0.95, 0.99)


@dataclass
class RunAnalysis:
    metadata: dict
    track_path: Path
    status: str
    message: str
    row_count: int = 0
    transition_count: int = 0
    used_count: int = 0
    downhill_count: int = 0
    uphill_count: int = 0
    observed_acceptance: float = math.nan
    mean_uphill_delta: float = math.nan
    median_uphill_delta: float = math.nan
    max_probability_error: float = math.nan
    initial_params_match: bool | None = None
    best_cost_error: float = math.nan
    deltas: list[float] | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate T0 per run and per lambda from a SQLite results database "
            "and its referenced SA tracking CSV files."
        )
    )
    parser.add_argument("--db", required=True, type=Path, help="SQLite database path")
    parser.add_argument(
        "--table", default="results", help="SQLite table name (default: results)"
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Root used for relative track_file_name paths (default: current directory)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("t0_analysis"),
        help="Directory for generated reports (default: t0_analysis)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=200,
        help="Initial transitions used per run; 0 uses all transitions (default: 200)",
    )
    parser.add_argument(
        "--targets",
        type=float,
        nargs="+",
        default=list(DEFAULT_TARGETS),
        help="Target initial acceptance rates (default: 0.80 0.85 0.90 0.95)",
    )
    parser.add_argument(
        "--primary-target",
        type=float,
        default=0.85,
        help="Target displayed in the text summary (default: 0.85)",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=5000,
        help="Bootstrap replicates for the mean per-run T0 (default: 5000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260903,
        help="Bootstrap random seed (default: 20260903)",
    )
    parser.add_argument(
        "--probability-tolerance",
        type=float,
        default=1e-10,
        help="Warning threshold for Metropolis probability verification",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit unsuccessfully when any tracking file is missing or invalid",
    )
    args = parser.parse_args()

    if args.window < 0:
        parser.error("--window must be nonnegative")
    if args.bootstrap < 0:
        parser.error("--bootstrap must be nonnegative")
    if not args.targets:
        parser.error("at least one target acceptance rate is required")
    if any(not 0.0 < target < 1.0 for target in args.targets):
        parser.error("all target acceptance rates must lie strictly between 0 and 1")
    if not 0.0 < args.primary_target < 1.0:
        parser.error("--primary-target must lie strictly between 0 and 1")

    args.targets = tuple(sorted(set(args.targets + [args.primary_target])))
    return args


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def parse_vector(value: object) -> list[float] | None:
    if value is None or str(value).strip() == "":
        return None
    parsed = ast.literal_eval(str(value))
    if not isinstance(parsed, (list, tuple)):
        raise ValueError("parameter vector is not a list or tuple")
    return [float(item) for item in parsed]


def vectors_close(
    first: Sequence[float] | None,
    second: Sequence[float] | None,
    tolerance: float = 1e-10,
) -> bool | None:
    if first is None or second is None:
        return None
    if len(first) != len(second):
        return False
    return all(
        math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)
        for a, b in zip(first, second)
    )


def resolve_track_path(track_name: str, db_path: Path, root: Path) -> Path:
    supplied = Path(track_name).expanduser()
    candidates = []
    if supplied.is_absolute():
        candidates.append(supplied)
    else:
        candidates.extend((root / supplied, db_path.parent / supplied))

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return candidates[0].resolve() if candidates else supplied.resolve()


def read_database(db_path: Path, table_name: str = "results") -> list[dict]:
    connection = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        table = connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (table_name,),
        ).fetchone()
        if table is None:
            raise ValueError(f"database has no {table_name!r} table")
        quoted_table = '"' + table_name.replace('"', '""') + '"'
        columns = {
            row[1] for row in connection.execute(f"PRAGMA table_info({quoted_table})")
        }
        required = {"run_id", "linear_penalization", "track_file_name"}
        missing = required - columns
        if missing:
            raise ValueError(f"{table_name!r} table is missing columns: {sorted(missing)}")
        return [dict(row) for row in connection.execute(f"SELECT * FROM {quoted_table} ORDER BY run_id")]
    finally:
        connection.close()


def metropolis_probability(delta: float, temperature: float) -> float:
    if delta <= 0.0:
        return 1.0
    if temperature <= 0.0:
        return 0.0
    exponent = -delta / temperature
    return 0.0 if exponent < -745.0 else math.exp(exponent)


def analyze_run(metadata: dict, db_path: Path, root: Path, window: int) -> RunAnalysis:
    track_name = str(metadata.get("track_file_name") or "").strip()
    track_path = resolve_track_path(track_name, db_path, root)
    result = RunAnalysis(
        metadata=metadata,
        track_path=track_path,
        status="error",
        message="",
    )
    if not track_name:
        result.message = "empty track_file_name"
        return result
    if not track_path.is_file():
        result.status = "missing"
        result.message = "tracking CSV not found"
        return result

    try:
        with track_path.open(newline="", encoding="utf-8-sig") as stream:
            rows = list(csv.DictReader(stream))
        required = {"temp", "cost", "accepted", "param"}
        if not rows:
            raise ValueError("tracking CSV has no data rows")
        missing = required - set(rows[0])
        if missing:
            raise ValueError(f"tracking CSV is missing columns: {sorted(missing)}")

        result.row_count = len(rows)
        first_params = parse_vector(rows[0].get("param"))
        database_params = parse_vector(metadata.get("initial_params"))
        result.initial_params_match = vectors_close(first_params, database_params)

        current_cost = float(rows[0]["cost"])
        deltas: list[float] = []
        accepted_flags: list[bool] = []
        probability_errors: list[float] = []

        for row in rows[1:]:
            proposed_cost = float(row["cost"])
            temperature = float(row["temp"])
            accepted = parse_bool(row["accepted"])
            delta = proposed_cost - current_cost
            deltas.append(delta)
            accepted_flags.append(accepted)

            recorded_probability = str(row.get("accept_prob", "")).strip()
            if recorded_probability:
                expected = metropolis_probability(delta, temperature)
                probability_errors.append(abs(float(recorded_probability) - expected))

            if accepted:
                current_cost = proposed_cost

        result.transition_count = len(deltas)
        used_deltas = deltas if window == 0 else deltas[:window]
        used_acceptance = accepted_flags if window == 0 else accepted_flags[:window]
        if not used_deltas:
            raise ValueError("no transitions available after the initial row")

        uphill = [delta for delta in used_deltas if delta > 0.0]
        result.used_count = len(used_deltas)
        result.downhill_count = len(used_deltas) - len(uphill)
        result.uphill_count = len(uphill)
        result.observed_acceptance = sum(used_acceptance) / len(used_acceptance)
        if uphill:
            result.mean_uphill_delta = statistics.fmean(uphill)
            result.median_uphill_delta = statistics.median(uphill)
        result.max_probability_error = max(probability_errors, default=math.nan)
        result.deltas = used_deltas

        database_best = metadata.get("best_cost")
        if database_best is not None:
            csv_best = min(float(row["cost"]) for row in rows)
            result.best_cost_error = abs(csv_best - float(database_best))

        messages = []
        if result.initial_params_match is False:
            messages.append("first CSV parameters do not match database initial_params")
        result.status = "warning" if messages else "ok"
        result.message = "; ".join(messages)
        return result
    except Exception as error:
        result.status = "error"
        result.message = f"{type(error).__name__}: {error}"
        return result


def empirical_acceptance(deltas: Sequence[float], temperature: float) -> float:
    if not deltas:
        return math.nan
    downhill = sum(delta <= 0.0 for delta in deltas)
    uphill_probability = sum(
        metropolis_probability(delta, temperature)
        for delta in deltas
        if delta > 0.0
    )
    return (downhill + uphill_probability) / len(deltas)


def empirical_t0(deltas: Sequence[float], target: float) -> float:
    if not deltas:
        return math.nan
    downhill_fraction = sum(delta <= 0.0 for delta in deltas) / len(deltas)
    uphill = [delta for delta in deltas if delta > 0.0]
    if not uphill or target <= downhill_fraction:
        return math.nan

    low = 0.0
    high = max(statistics.fmean(uphill), 1e-12)
    while empirical_acceptance(deltas, high) < target:
        high *= 2.0
        if not math.isfinite(high):
            return math.nan

    for _ in range(80):
        midpoint = (low + high) / 2.0
        if empirical_acceptance(deltas, midpoint) < target:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2.0


def paper_t0(deltas: Sequence[float], target: float) -> float:
    if not deltas:
        return math.nan
    uphill = [delta for delta in deltas if delta > 0.0]
    m2 = len(uphill)
    m1 = len(deltas) - m2
    if not uphill:
        return math.nan
    denominator = m2 * target - m1 * (1.0 - target)
    if denominator <= 0.0 or denominator >= m2:
        return math.nan
    return statistics.fmean(uphill) / math.log(m2 / denominator)


def percentile(values: Sequence[float], probability: float) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    if len(clean) == 1:
        return clean[0]
    position = probability * (len(clean) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    fraction = position - lower
    return clean[lower] * (1.0 - fraction) + clean[upper] * fraction


def bootstrap_mean_interval(
    values: Sequence[float], replicates: int, seed: int
) -> tuple[float, float]:
    clean = [value for value in values if math.isfinite(value)]
    if not clean or replicates <= 0:
        return math.nan, math.nan
    if len(clean) == 1:
        return clean[0], clean[0]
    generator = random.Random(seed)
    means = [
        statistics.fmean(generator.choices(clean, k=len(clean)))
        for _ in range(replicates)
    ]
    return percentile(means, 0.025), percentile(means, 0.975)


def float_text(value: object, digits: int = 12) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return "" if not math.isfinite(number) else f"{number:.{digits}g}"


def write_csv(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: float_text(row.get(field)) for field in fields})


def run_diagnostic_row(analysis: RunAnalysis, probability_tolerance: float) -> dict:
    metadata = analysis.metadata
    probability_valid = (
        math.isfinite(analysis.max_probability_error)
        and analysis.max_probability_error <= probability_tolerance
    )
    return {
        "run_id": metadata.get("run_id"),
        "lambda": metadata.get("linear_penalization"),
        "track_file_name": metadata.get("track_file_name"),
        "resolved_track_path": str(analysis.track_path),
        "status": analysis.status,
        "message": analysis.message,
        "row_count": analysis.row_count,
        "transition_count": analysis.transition_count,
        "used_count": analysis.used_count,
        "downhill_count": analysis.downhill_count,
        "uphill_count": analysis.uphill_count,
        "observed_acceptance": analysis.observed_acceptance,
        "mean_uphill_delta": analysis.mean_uphill_delta,
        "median_uphill_delta": analysis.median_uphill_delta,
        "max_probability_error": analysis.max_probability_error,
        "probability_valid": probability_valid,
        "initial_params_match": analysis.initial_params_match,
        "best_cost_error": analysis.best_cost_error,
        "initial_temp": metadata.get("initial_temp"),
        "min_temp": metadata.get("min_temp"),
        "cooling_rate": metadata.get("cooling_rate"),
        "perturbation": metadata.get("perturbation"),
        "gen_seed": metadata.get("gen_seed"),
    }


def write_reports(args: argparse.Namespace, analyses: list[RunAnalysis]) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)

    diagnostic_fields = (
        "run_id",
        "lambda",
        "track_file_name",
        "resolved_track_path",
        "status",
        "message",
        "row_count",
        "transition_count",
        "used_count",
        "downhill_count",
        "uphill_count",
        "observed_acceptance",
        "mean_uphill_delta",
        "median_uphill_delta",
        "max_probability_error",
        "probability_valid",
        "initial_params_match",
        "best_cost_error",
        "initial_temp",
        "min_temp",
        "cooling_rate",
        "perturbation",
        "gen_seed",
    )
    diagnostic_rows = [
        run_diagnostic_row(analysis, args.probability_tolerance)
        for analysis in analyses
    ]
    write_csv(args.output_dir / "t0_diagnostics.csv", diagnostic_rows, diagnostic_fields)

    per_run_fields = (
        "run_id",
        "lambda",
        "target_acceptance",
        "empirical_t0",
        "paper_t0",
        "used_count",
        "downhill_fraction",
        "mean_uphill_delta",
        "observed_acceptance",
        "status",
        "track_file_name",
    )
    per_run_rows = []
    estimates_by_group: dict[tuple[float, float], list[float]] = defaultdict(list)
    deltas_by_lambda: dict[float, list[float]] = defaultdict(list)

    for analysis in analyses:
        if analysis.status not in {"ok", "warning"} or not analysis.deltas:
            continue
        lambda_value = float(analysis.metadata["linear_penalization"])
        deltas_by_lambda[lambda_value].extend(analysis.deltas)
        for target in args.targets:
            estimate = empirical_t0(analysis.deltas, target)
            approximation = paper_t0(analysis.deltas, target)
            if math.isfinite(estimate):
                estimates_by_group[(lambda_value, target)].append(estimate)
            per_run_rows.append(
                {
                    "run_id": analysis.metadata.get("run_id"),
                    "lambda": lambda_value,
                    "target_acceptance": target,
                    "empirical_t0": estimate,
                    "paper_t0": approximation,
                    "used_count": analysis.used_count,
                    "downhill_fraction": analysis.downhill_count / analysis.used_count,
                    "mean_uphill_delta": analysis.mean_uphill_delta,
                    "observed_acceptance": analysis.observed_acceptance,
                    "status": analysis.status,
                    "track_file_name": analysis.metadata.get("track_file_name"),
                }
            )
    write_csv(args.output_dir / "t0_per_run.csv", per_run_rows, per_run_fields)

    lambda_fields = (
        "lambda",
        "target_acceptance",
        "valid_run_count",
        "pooled_transition_count",
        "pooled_empirical_t0",
        "pooled_paper_t0",
        "mean_run_t0",
        "median_run_t0",
        "sd_run_t0",
        "min_run_t0",
        "p10_run_t0",
        "p90_run_t0",
        "max_run_t0",
        "bootstrap_mean_ci_low",
        "bootstrap_mean_ci_high",
    )
    lambda_rows = []
    for lambda_value in sorted(deltas_by_lambda):
        pooled_deltas = deltas_by_lambda[lambda_value]
        for target in args.targets:
            estimates = estimates_by_group[(lambda_value, target)]
            group_seed = args.seed + int(round(lambda_value * 1000)) + int(round(target * 100))
            ci_low, ci_high = bootstrap_mean_interval(
                estimates, args.bootstrap, group_seed
            )
            lambda_rows.append(
                {
                    "lambda": lambda_value,
                    "target_acceptance": target,
                    "valid_run_count": len(estimates),
                    "pooled_transition_count": len(pooled_deltas),
                    "pooled_empirical_t0": empirical_t0(pooled_deltas, target),
                    "pooled_paper_t0": paper_t0(pooled_deltas, target),
                    "mean_run_t0": statistics.fmean(estimates) if estimates else math.nan,
                    "median_run_t0": statistics.median(estimates) if estimates else math.nan,
                    "sd_run_t0": statistics.stdev(estimates) if len(estimates) > 1 else math.nan,
                    "min_run_t0": min(estimates) if estimates else math.nan,
                    "p10_run_t0": percentile(estimates, 0.10),
                    "p90_run_t0": percentile(estimates, 0.90),
                    "max_run_t0": max(estimates) if estimates else math.nan,
                    "bootstrap_mean_ci_low": ci_low,
                    "bootstrap_mean_ci_high": ci_high,
                }
            )
    write_csv(args.output_dir / "t0_by_lambda.csv", lambda_rows, lambda_fields)

    ok_count = sum(analysis.status in {"ok", "warning"} for analysis in analyses)
    missing_count = sum(analysis.status == "missing" for analysis in analyses)
    error_count = sum(analysis.status == "error" for analysis in analyses)
    warning_count = sum(analysis.status == "warning" for analysis in analyses)
    primary_rows = [
        row
        for row in lambda_rows
        if math.isclose(
            float(row["target_acceptance"]), args.primary_target, abs_tol=1e-12
        )
    ]
    summary_lines = [
        "Simulated-annealing initial-temperature analysis",
        "================================================",
        f"Script version: {SCRIPT_VERSION}",
        f"Database: {args.db.resolve()}",
        f"Table: {args.table}",
        f"Database SHA-256: {sha256_file(args.db)}",
        f"Track root: {args.root.resolve()}",
        f"Initial transition window: {'all' if args.window == 0 else args.window}",
        f"Target acceptance rates: {', '.join(f'{target:.2f}' for target in args.targets)}",
        f"Bootstrap replicates: {args.bootstrap}",
        f"Bootstrap seed: {args.seed}",
        "",
        f"Database runs: {len(analyses)}",
        f"Usable runs: {ok_count}",
        f"Warnings: {warning_count}",
        f"Missing tracking files: {missing_count}",
        f"Invalid runs: {error_count}",
        "",
        f"Results for target acceptance {args.primary_target:.2f}",
        "lambda,runs,pooled_t0,mean_run_t0,sd_run_t0,p10_run_t0,p90_run_t0,mean_ci_low,mean_ci_high",
    ]
    for row in primary_rows:
        summary_lines.append(
            ",".join(
                (
                    float_text(row["lambda"]),
                    float_text(row["valid_run_count"]),
                    float_text(row["pooled_empirical_t0"]),
                    float_text(row["mean_run_t0"]),
                    float_text(row["sd_run_t0"]),
                    float_text(row["p10_run_t0"]),
                    float_text(row["p90_run_t0"]),
                    float_text(row["bootstrap_mean_ci_low"]),
                    float_text(row["bootstrap_mean_ci_high"]),
                )
            )
        )
    summary_lines.extend(
        (
            "",
            "Notes",
            "-----",
            "The first CSV row is treated as the initial state, not as a proposed transition.",
            "Empirical T0 solves the acceptance equation using every positive cost difference.",
            "Paper T0 replaces the positive cost differences by their arithmetic mean.",
            "Bootstrap intervals describe uncertainty in the mean of the per-run T0 estimates.",
            "Use t0_diagnostics.csv to inspect missing files and consistency checks.",
        )
    )
    (args.output_dir / "t0_summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )

    manifest = {
        "script_version": SCRIPT_VERSION,
        "database": str(args.db.resolve()),
        "table": args.table,
        "database_sha256": sha256_file(args.db),
        "track_root": str(args.root.resolve()),
        "window": args.window,
        "targets": args.targets,
        "primary_target": args.primary_target,
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_seed": args.seed,
        "probability_tolerance": args.probability_tolerance,
        "database_run_count": len(analyses),
        "usable_run_count": ok_count,
        "warning_count": warning_count,
        "missing_file_count": missing_count,
        "error_count": error_count,
    }
    (args.output_dir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    args = parse_args()
    if not args.db.is_file():
        print(f"error: database not found: {args.db}", file=sys.stderr)
        return 2

    try:
        metadata_rows = read_database(args.db, args.table)
    except Exception as error:
        print(f"error: {type(error).__name__}: {error}", file=sys.stderr)
        return 2

    analyses = [
        analyze_run(metadata, args.db, args.root.resolve(), args.window)
        for metadata in metadata_rows
    ]
    write_reports(args, analyses)

    missing_or_invalid = sum(
        analysis.status in {"missing", "error"} for analysis in analyses
    )
    print(f"Analyzed {len(analyses) - missing_or_invalid}/{len(analyses)} runs.")
    print(f"Reports written to {args.output_dir.resolve()}")
    if missing_or_invalid:
        print(
            f"{missing_or_invalid} runs were missing or invalid; see t0_diagnostics.csv.",
            file=sys.stderr,
        )
    return 1 if args.strict and missing_or_invalid else 0


if __name__ == "__main__":
    raise SystemExit(main())
