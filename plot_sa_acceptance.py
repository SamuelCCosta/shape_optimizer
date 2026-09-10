#!/usr/bin/env python3
"""Plot moving SA acceptance rates or probabilities against temperature.

Examples:
    python plot_sa_acceptance.py --db experiments.db --table low_temp_n1
    python plot_sa_acceptance.py --db experiments.db --run-ids 1 2 --window 100
    python plot_sa_acceptance.py --db experiments.db --table low_temp_n1 \
        --metric uphill-acceptance
    python plot_sa_acceptance.py --db experiments.db --table low_temp_n1 \
        --metric uphill-probability
    python plot_sa_acceptance.py --csv track_csv/experiments/low_temp_n1/1.csv \
        --window 200 --output plots/sa_acceptance.png

The first CSV data row is the initial state, not a move. Each plotted point
uses exactly the last ``window`` moves, including the current move, and is
plotted at the current move's temperature. Partial windows are omitted to
avoid noisy estimates at the high-temperature end of the curve.
Windows follow CSV order and never cross run boundaries. The default metric
uses accepted/all moves. Both uphill metrics first filter to ``accept_prob != 1``:
uphill-acceptance measures accepted/worsening moves; uphill-probability averages
the CSV accept_prob values. Their windows count only worsening moves, and their
y axes display percentages. Requires Matplotlib to render the figure.
"""

from __future__ import annotations

import argparse
import csv
import math
import sqlite3
from collections import deque
from pathlib import Path

from estimate_sa_initial_temperature import (
    parse_bool,
    read_database,
    resolve_track_path,
)


METRICS = ("acceptance", "uphill-acceptance", "uphill-probability")


def moving_acceptance(
    track_path: Path, window: int = 200, metric: str = "acceptance"
) -> tuple[list[float], list[float]]:
    """Return temperatures and full-window means, expressed as fractions (0–1).

    Uphill modes filter to accept_prob != 1 before counting the window.
    The first result is at eligible proposal ``window``; partial windows are
    omitted. Too few eligible proposals raises ValueError.
    Nonpositive or nonfinite proposal temperatures are rejected because they
    cannot be represented on a logarithmic temperature axis.
    """
    if window < 1:
        raise ValueError("window must be positive")
    if metric not in METRICS:
        raise ValueError(f"unknown metric: {metric!r}")
    uphill_only = metric != "acceptance"
    temperatures: list[float] = []
    rates: list[float] = []
    recent: deque[float] = deque()
    window_sum = 0.0
    move_count = 0
    with track_path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {"temp"}
        if metric != "uphill-probability":
            required.add("accepted")
        if uphill_only:
            required.add("accept_prob")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{track_path}: missing columns: {sorted(missing)}")
        next(reader, None)  # Initial state, logged as accepted but not a proposal.
        for row in reader:
            try:
                probability = 1.0
                if uphill_only:
                    probability = float(row["accept_prob"])
                    if not math.isfinite(probability) or not 0 <= probability <= 1:
                        raise ValueError("accept_prob must be finite and between 0 and 1")
                    if probability == 1:
                        continue
                temperature = float(row["temp"])
                if not math.isfinite(temperature) or temperature <= 0:
                    raise ValueError("temperature must be finite and positive")
                value = probability if metric == "uphill-probability" else float(parse_bool(row["accepted"]))
            except (TypeError, ValueError) as error:
                raise ValueError(f"{track_path}: line {reader.line_num}: {error}") from error
            if len(recent) == window:
                window_sum -= recent.popleft()
            recent.append(value)
            window_sum += value
            move_count += 1
            if len(recent) == window:
                temperatures.append(temperature)
                # Bound tiny floating-point roundoff when averaging probabilities.
                rates.append(max(0.0, min(1.0, window_sum / window)))
    if not temperatures:
        raise ValueError(
            f"{track_path}: only {move_count} {'worsening ' if uphill_only else ''}moves after the initial state; "
            f"a full window requires {window}. Choose a smaller --window."
        )
    return temperatures, rates


def plot_acceptance(
    tracks: list[tuple[str, Path]], window: int = 200, metric: str = "acceptance"
):
    """Create and return a Matplotlib figure with one curve per tracking CSV."""
    import matplotlib.pyplot as plt

    if not tracks:
        raise ValueError("no tracking files selected")
    curves = [(label, *moving_acceptance(path, window, metric)) for label, path in tracks]
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    for label, temperatures, rates in curves:
        ax.plot(temperatures, rates, linewidth=1.3, alpha=0.85, label=label)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("Temperature (log scale)")
    labels = {
        "acceptance": ("Acceptance rate (accepted / all moves)", "SA moving acceptance rate", "moves"),
        "uphill-acceptance": ("Worsening moves accepted (%)", "Observed acceptance of worsening moves", "worsening moves"),
        "uphill-probability": ("Mean acceptance probability (%)", "Mean CSV acceptance probability of worsening moves", "worsening moves"),
    }
    ylabel, title, window_unit = labels[metric]
    ax.set_ylabel(ylabel)
    if metric != "acceptance":
        from matplotlib.ticker import PercentFormatter
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_ylim(0, 1)
    ax.set_title(f"{title}\nFull trailing windows of {window} {window_unit}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize="small")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--db", type=Path, help="SQLite database, as in estimate_sa_initial_temperature.py")
    source.add_argument("--csv", type=Path, nargs="+", help="Tracking CSV files to plot directly")
    parser.add_argument("--table", default="results", help="Database table (default: results)")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Root for relative database tracking paths")
    parser.add_argument("--run-ids", type=int, nargs="+", help="Database runs to plot (default: all)")
    parser.add_argument("--metric", choices=METRICS, default="acceptance", help="Observed acceptance of all moves (default), observed acceptance of worsening moves, or mean CSV probability of worsening moves")
    parser.add_argument("--window", type=int, default=200, help="Full trailing window in eligible moves (only worsening moves for uphill metrics); partial windows omitted (default: 200)")
    parser.add_argument("--output", type=Path, help="Save figure instead of displaying it; e.g. plots/sa_acceptance.png")
    args = parser.parse_args()
    if args.window < 1:
        parser.error("--window must be positive")
    if args.csv and args.run_ids:
        parser.error("--run-ids requires --db")

    try:
        if args.csv:
            tracks = [(str(path), path) for path in args.csv]
        else:
            rows = read_database(args.db, args.table)
            if args.run_ids:
                missing = set(args.run_ids) - {row["run_id"] for row in rows}
                if missing:
                    raise ValueError(f"unknown run IDs: {sorted(missing)}")
                rows = [row for row in rows if row["run_id"] in args.run_ids]
            tracks = []
            for row in rows:
                track_name = str(row.get("track_file_name") or "").strip()
                if not track_name:
                    raise ValueError(f"run {row['run_id']}: empty track_file_name")
                label = f"Run {row['run_id']} (λ={row['linear_penalization']:g})"
                tracks.append((label, resolve_track_path(track_name, args.db, args.root)))

        if args.output:
            import matplotlib
            matplotlib.use("Agg")
        fig = plot_acceptance(tracks, args.window, args.metric)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(args.output, dpi=180)
            print(f"Saved {args.output} ({len(tracks)} runs, window={args.window}, metric={args.metric})")
        else:
            import matplotlib.pyplot as plt
            plt.show()
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
