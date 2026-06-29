from __future__ import annotations

import csv
import gzip
import math
import re
import sys
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# Your data folder. The r"..." raw string is important for Windows paths.
DATA_FOLDER = Path(r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5 ToF IV R\G5")

# Your program can be placed anywhere, e.g.:
# C:\Users\jsjan\Downloads\SSC Data&Codes\livetime.py
# The input data are taken from DATA_FOLDER above.

# To process specific files only, put their paths here. If this list is empty,
# all .zip, .csv, and .csv.gz files in DATA_FOLDER are processed.
INPUT_FILES: list[Path] = []
# Example:
# INPUT_FILES = [
#     Path(r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5 ToF IV R\G5\lifetime_50avg.zip"),
#     Path(r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5 ToF IV R\G5\lifetime_200avg.zip"),
# ]

# Where to save figures and the table. If None, a folder named
# "lifetime_results" is created inside DATA_FOLDER or beside the selected files.
OUTPUT_FOLDER: Optional[Path] = None

# Read and analyse this time window around the first decay.
# Your uploaded files show the relevant decay before the next 1 kHz cycle, so
# 0...500 us is a good default. Change ANALYSIS_END_US if your trace differs.
PRETRIGGER_US = 20.0
ANALYSIS_END_US = 500.0

# Baseline is estimated from the final part of the analysis window.
# With ANALYSIS_END_US = 500 and BASELINE_WINDOW_US = 50, the baseline is the
# median signal from 450...500 us.
BASELINE_WINDOW_US = 50.0

# Average this many raw oscilloscope samples into one plotted/fitted point.
# The raw files contain 25 million points; binning makes the plot fast and
# greatly reduces noise without changing the microsecond-scale lifetime.
BIN_POINTS = 1000

# Automatic decay-start detection. The code first detects the high transient
# plateau, then starts the fit when the signal has fallen below this fraction
# of its peak. For a trace with no plateau, this will be close to the peak.
DECAY_START_FRACTION = 0.98

# Initial-slope fitting.
# Both methods use the same immediate points after the detected decay start,
# instead of selecting a time interval or signal-fraction interval.
# The code starts with INITIAL_SLOPE_MIN_POINTS and, only if noise makes the
# initial fitted slope non-negative, extends up to INITIAL_SLOPE_MAX_POINTS.
INITIAL_SLOPE_MIN_POINTS = 8
INITIAL_SLOPE_MAX_POINTS = 20

# Skip this many binned points after the detected start. Usually 0.
# Use 1 or 2 only if the very first point is distorted by trigger/ringing.
INITIAL_SLOPE_SKIP_POINTS = 3

# Manual overrides. Leave as None for automatic detection. Values are in raw
# oscilloscope time in microseconds, not in shifted plot time.
MANUAL_DECAY_START_US: Optional[float] = None
MANUAL_BASELINE_V: Optional[float] = None

# Show the figures interactively after saving them. Useful when double-clicking.
SHOW_PLOTS = True

# Keep the console window open at the end when double-clicking in Windows.
PAUSE_AT_END = True


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclass
class LifetimeResult:
    file: str
    baseline_v: float
    polarity: int
    raw_decay_start_us: float
    peak_signal_mV: float
    tau_method1_us: float
    tau_method2_us: float
    method1_slope_mV_per_us: float
    method2_slope_per_us: float
    method1_points: int
    method2_points: int
    figure_path: str


# ---------------------------------------------------------------------------
# FILE HANDLING
# ---------------------------------------------------------------------------

def _choose_csv_inside_zip(zip_path: Path) -> str:
    """Return the first CSV-like member in a zip file."""
    with zipfile.ZipFile(zip_path) as zf:
        names = [n for n in zf.namelist() if not n.endswith("/")]
        csv_names = [n for n in names if n.lower().endswith((".csv", ".txt"))]
        if not csv_names:
            raise ValueError(f"No CSV/TXT file found inside {zip_path}")
        return csv_names[0]


@contextmanager
def open_text_or_binary_csv(path: Path) -> Iterator[Iterable[bytes]]:
    """Open a CSV, CSV.GZ, or ZIP containing a CSV as a binary file-like object."""
    suffixes = "".join(path.suffixes).lower()
    if path.suffix.lower() == ".zip":
        zf = zipfile.ZipFile(path)
        member = _choose_csv_inside_zip(path)
        try:
            with zf.open(member, "r") as f:
                yield f
        finally:
            zf.close()
    elif suffixes.endswith(".csv.gz") or path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as f:
            yield f
    else:
        with open(path, "rb") as f:
            yield f


def find_waveform_header_row(path: Path, max_lines: int = 200) -> int:
    """Find the row containing the waveform column names, normally TIME,CH1."""
    with open_text_or_binary_csv(path) as f:
        for row_number, raw_line in enumerate(f):
            if row_number > max_lines:
                break
            line = raw_line.decode("utf-8", errors="replace").strip().replace("\ufeff", "")
            parts = [p.strip().strip('"').upper() for p in line.split(",")]
            if len(parts) >= 2 and parts[0] == "TIME":
                return row_number
    raise ValueError(
        f"Could not find a TIME column in the first {max_lines} lines of {path}."
    )


def read_waveform_window(path: Path, t_min_s: float, t_max_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Read only the requested time window from a large Tektronix CSV file."""
    header_row = find_waveform_header_row(path)
    chunks: list[pd.DataFrame] = []

    with open_text_or_binary_csv(path) as f:
        reader = pd.read_csv(
            f,
            skiprows=header_row,
            usecols=[0, 1],
            chunksize=1_000_000,
            low_memory=False,
        )

        for chunk in reader:
            # Normalize column names to TIME and CH1 even if the channel is named differently.
            chunk = chunk.rename(columns={chunk.columns[0]: "TIME", chunk.columns[1]: "CH1"})
            chunk["TIME"] = pd.to_numeric(chunk["TIME"], errors="coerce")
            chunk["CH1"] = pd.to_numeric(chunk["CH1"], errors="coerce")
            chunk = chunk.dropna(subset=["TIME", "CH1"])
            if chunk.empty:
                continue

            in_window = (chunk["TIME"] >= t_min_s) & (chunk["TIME"] <= t_max_s)
            if in_window.any():
                chunks.append(chunk.loc[in_window, ["TIME", "CH1"]])

            # Tektronix files are time sorted. Stop once we have passed t_max_s.
            if float(chunk["TIME"].iloc[-1]) > t_max_s:
                break

    if not chunks:
        raise ValueError(f"No data found between {t_min_s:g} s and {t_max_s:g} s in {path}")

    df = pd.concat(chunks, ignore_index=True)
    return df["TIME"].to_numpy(dtype=float), df["CH1"].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# ANALYSIS
# ---------------------------------------------------------------------------

def select_initial_slope_region(
    t_rel_s: np.ndarray,
    signal_v: np.ndarray,
    min_points: int,
    max_points: int,
    skip_points: int = 0,
) -> np.ndarray:
    """
    Select the immediate first valid points after the detected decay start.

    This makes the two lifetime estimates local initial-slope estimates:
    - Method 1 fits signal vs time over these initial points.
    - Method 2 fits log(signal) vs time over the same initial points.

    The selection is based on consecutive point count, not a time interval and
    not high/low signal-fraction thresholds. To avoid a single noisy point
    producing a positive initial slope, the selector uses the shortest prefix
    between min_points and max_points whose linear and log slopes are both
    negative.
    """
    valid = np.where(
        (t_rel_s >= 0.0)
        & np.isfinite(t_rel_s)
        & np.isfinite(signal_v)
        & (signal_v > 0.0)
    )[0]

    if len(valid) <= skip_points:
        return np.zeros_like(t_rel_s, dtype=bool)

    valid = valid[skip_points:]
    max_n = min(max_points, len(valid))
    min_n = min(min_points, max_n)

    chosen = valid[:min_n]
    for n in range(min_n, max_n + 1):
        candidate = valid[:n]
        linear_slope = np.polyfit(t_rel_s[candidate], signal_v[candidate], 1)[0]
        log_slope = np.polyfit(t_rel_s[candidate], np.log(signal_v[candidate]), 1)[0]
        chosen = candidate
        if linear_slope < 0.0 and log_slope < 0.0:
            break

    mask = np.zeros_like(t_rel_s, dtype=bool)
    mask[chosen] = True
    return mask

def bin_average(t_s: np.ndarray, y: np.ndarray, points_per_bin: int) -> tuple[np.ndarray, np.ndarray]:
    """Average adjacent points so the 25M-point trace can be plotted and fitted quickly."""
    if points_per_bin <= 1 or len(t_s) < 2 * points_per_bin:
        return t_s, y
    n_bins = len(t_s) // points_per_bin
    trimmed = n_bins * points_per_bin
    t_binned = t_s[:trimmed].reshape(n_bins, points_per_bin).mean(axis=1)
    y_binned = y[:trimmed].reshape(n_bins, points_per_bin).mean(axis=1)
    return t_binned, y_binned


def estimate_baseline(t_s: np.ndarray, voltage_v: np.ndarray) -> float:
    if MANUAL_BASELINE_V is not None:
        return float(MANUAL_BASELINE_V)

    start_s = (ANALYSIS_END_US - BASELINE_WINDOW_US) * 1e-6
    end_s = ANALYSIS_END_US * 1e-6
    mask = (t_s >= start_s) & (t_s <= end_s)
    if mask.sum() >= 10:
        return float(np.median(voltage_v[mask]))

    # Fallback: median of final 10% of the loaded window.
    n_tail = max(10, len(voltage_v) // 10)
    return float(np.median(voltage_v[-n_tail:]))


def determine_polarity(t_s: np.ndarray, voltage_v: np.ndarray, baseline_v: float) -> int:
    """Return +1 for positive-going transient, -1 for negative-going transient."""
    mask = (t_s >= 0.0) & (t_s <= min(50e-6, ANALYSIS_END_US * 1e-6))
    if mask.sum() < 10:
        mask = np.ones_like(t_s, dtype=bool)
    deviation = voltage_v[mask] - baseline_v
    idx = int(np.nanargmax(np.abs(deviation)))
    sign = int(np.sign(deviation[idx]))
    return sign if sign != 0 else 1


def detect_decay_start(t_s: np.ndarray, signal_v: np.ndarray) -> float:
    """Detect where the actual signal fall starts and return raw scope time in seconds."""
    if MANUAL_DECAY_START_US is not None:
        return float(MANUAL_DECAY_START_US) * 1e-6

    mask = (t_s >= 0.0) & (t_s <= min(50e-6, ANALYSIS_END_US * 1e-6))
    if mask.sum() < 10:
        return 0.0

    peak_v = float(np.percentile(signal_v[mask], 99.5))
    if not np.isfinite(peak_v) or peak_v <= 0:
        return 0.0

    # Find first approach to the high transient region.
    high_indices = np.where((t_s >= 0.0) & (signal_v > 0.80 * peak_v))[0]
    if len(high_indices) == 0:
        return 0.0
    first_high = high_indices[0]

    # Start the lifetime fit when the high plateau has begun to fall.
    below_after_high = np.where(
        (np.arange(len(t_s)) > first_high)
        & (t_s >= 0.0)
        & (signal_v < DECAY_START_FRACTION * peak_v)
    )[0]
    if len(below_after_high) == 0:
        return float(t_s[first_high])

    return float(t_s[below_after_high[0]])


def signal_at_decay_start(t_rel_s: np.ndarray, signal_v: np.ndarray) -> float:
    """Estimate signal amplitude at t_rel = 0 for fit-window fractions."""
    mask = (t_rel_s >= 0.0) & (t_rel_s <= 1e-6) & (signal_v > 0.0)
    if mask.sum() >= 3:
        return float(np.median(signal_v[mask]))
    mask = (t_rel_s >= -2e-6) & (t_rel_s <= 2e-6) & (signal_v > 0.0)
    if mask.sum() >= 3:
        return float(np.percentile(signal_v[mask], 95))
    return float(np.nanmax(signal_v))


def select_fraction_fit_region(
    t_rel_s: np.ndarray,
    signal_v: np.ndarray,
    peak_v: float,
    high_fraction: float,
    low_fraction: float,
) -> np.ndarray:
    """Select the first monotonic-looking decay region between two signal fractions."""
    high = high_fraction * peak_v
    low = low_fraction * peak_v

    post = np.where((t_rel_s >= 0.0) & (signal_v > 0.0))[0]
    if len(post) == 0:
        return np.zeros_like(t_rel_s, dtype=bool)

    # Use only the first passage down to the low threshold. This avoids later noise.
    below_low = post[signal_v[post] < low]
    stop_index = below_low[0] if len(below_low) else post[-1]

    mask = (
        (t_rel_s >= 0.0)
        & (np.arange(len(t_rel_s)) <= stop_index)
        & (signal_v <= high)
        & (signal_v >= low)
    )
    return mask


def linear_method_1(t_rel_s: np.ndarray, signal_v: np.ndarray, mask: np.ndarray) -> tuple[float, float, float]:
    """Return tau_s, slope_v_per_s, intercept_v from the initial linear tangent."""
    if mask.sum() < 3:
        return math.nan, math.nan, math.nan
    slope, intercept = np.polyfit(t_rel_s[mask], signal_v[mask], 1)
    # For a local tangent at t = 0, tau is where the tangent reaches zero.
    tau_s = -intercept / slope if slope < 0 and intercept > 0 else math.nan
    return float(tau_s), float(slope), float(intercept)


def log_method_2(t_rel_s: np.ndarray, signal_v: np.ndarray, mask: np.ndarray) -> tuple[float, float, float]:
    """Return tau_s, slope_per_s, intercept_log_v from the initial log-slope."""
    fit = mask & (signal_v > 0.0)
    if fit.sum() < 3:
        return math.nan, math.nan, math.nan
    slope, intercept = np.polyfit(t_rel_s[fit], np.log(signal_v[fit]), 1)
    tau_s = -1.0 / slope if slope < 0 else math.nan
    return float(tau_s), float(slope), float(intercept)


def safe_stem(path: Path) -> str:
    stem = path.stem
    if stem.lower() == "csv" and len(path.suffixes) > 1:
        stem = Path(path.name.removesuffix(".gz")).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)


def analyse_one_file(path: Path, output_folder: Path) -> LifetimeResult:
    print(f"\nReading {path.name} ...")
    t_raw_s, v_raw = read_waveform_window(
        path,
        t_min_s=-PRETRIGGER_US * 1e-6,
        t_max_s=ANALYSIS_END_US * 1e-6,
    )

    baseline_v = estimate_baseline(t_raw_s, v_raw)
    polarity = determine_polarity(t_raw_s, v_raw, baseline_v)

    # Shift long-time baseline to 0 and flip negative-going traces upward.
    signal_raw_v = polarity * (v_raw - baseline_v)

    t_s, signal_v = bin_average(t_raw_s, signal_raw_v, BIN_POINTS)
    decay_start_s = detect_decay_start(t_s, signal_v)
    t_rel_s = t_s - decay_start_s
    peak_v = signal_at_decay_start(t_rel_s, signal_v)

    # Both fits are now based on the same immediate initial slope region.
    # This avoids fitting over an arbitrary time interval or amplitude interval.
    initial_slope_mask = select_initial_slope_region(
        t_rel_s,
        signal_v,
        INITIAL_SLOPE_MIN_POINTS,
        INITIAL_SLOPE_MAX_POINTS,
        INITIAL_SLOPE_SKIP_POINTS,
    )
    fit1_mask = initial_slope_mask
    fit2_mask = initial_slope_mask

    tau1_s, slope1_v_s, intercept1_v = linear_method_1(t_rel_s, signal_v, fit1_mask)
    tau2_s, slope2_per_s, intercept2_log = log_method_2(t_rel_s, signal_v, fit2_mask)

    figure_path = output_folder / f"{safe_stem(path)}_Figure76_lifetime.png"
    make_figure_76(
        path=path,
        figure_path=figure_path,
        t_rel_s=t_rel_s,
        signal_v=signal_v,
        fit1_mask=fit1_mask,
        fit2_mask=fit2_mask,
        tau1_s=tau1_s,
        slope1_v_s=slope1_v_s,
        intercept1_v=intercept1_v,
        tau2_s=tau2_s,
        slope2_per_s=slope2_per_s,
        intercept2_log=intercept2_log,
    )

    result = LifetimeResult(
        file=path.name,
        baseline_v=baseline_v,
        polarity=polarity,
        raw_decay_start_us=decay_start_s * 1e6,
        peak_signal_mV=peak_v * 1e3,
        tau_method1_us=tau1_s * 1e6 if np.isfinite(tau1_s) else math.nan,
        tau_method2_us=tau2_s * 1e6 if np.isfinite(tau2_s) else math.nan,
        method1_slope_mV_per_us=slope1_v_s * 1e-3 if np.isfinite(slope1_v_s) else math.nan,
        method2_slope_per_us=slope2_per_s * 1e-6 if np.isfinite(slope2_per_s) else math.nan,
        method1_points=int(fit1_mask.sum()),
        method2_points=int(fit2_mask.sum()),
        figure_path=str(figure_path),
    )

    print(f"  baseline        = {result.baseline_v:+.6g} V")
    print(f"  polarity        = {result.polarity:+d}")
    print(f"  decay start     = {result.raw_decay_start_us:.3f} us in raw scope time")
    print(f"  Method 1 tau    = {result.tau_method1_us:.3f} us")
    print(f"  Method 2 tau    = {result.tau_method2_us:.3f} us")
    print(f"  saved figure    = {figure_path}")
    return result


# ---------------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------------

def make_figure_76(
    path: Path,
    figure_path: Path,
    t_rel_s: np.ndarray,
    signal_v: np.ndarray,
    fit1_mask: np.ndarray,
    fit2_mask: np.ndarray,
    tau1_s: float,
    slope1_v_s: float,
    intercept1_v: float,
    tau2_s: float,
    slope2_per_s: float,
    intercept2_log: float,
) -> None:
    x_us = t_rel_s * 1e6
    y_mV = signal_v * 1e3

    # Plot only the important part of the decay. Keep positive values for log plot.
    plot_mask = (x_us >= 0.0) & (x_us <= min(ANALYSIS_END_US, 200.0))
    log_mask = plot_mask & (signal_v > 0.0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)

    # Figure 76a: linear-linear scale, straight-line intercept with baseline.
    ax = axes[0]
    ax.plot(x_us[plot_mask], y_mV[plot_mask], ".", markersize=2, label="baseline-shifted data")
    ax.axhline(0.0, linestyle="--", linewidth=1, label="baseline")
    ax.axvline(0.0, linestyle=":", linewidth=1, label="detected decay start")

    if np.isfinite(tau1_s) and np.isfinite(slope1_v_s) and np.isfinite(intercept1_v):
        fit_x_s = np.linspace(max(0.0, np.nanmin(t_rel_s[fit1_mask])), tau1_s * 1.08, 250)
        fit_y_mV = (slope1_v_s * fit_x_s + intercept1_v) * 1e3
        ax.plot(fit_x_s * 1e6, fit_y_mV, linewidth=2, label=f"initial tangent, tau = {tau1_s*1e6:.2f} $\mu$s")
        ax.axvline(tau1_s * 1e6, linestyle="--", linewidth=1)

    ax.set_xlabel("Time after decay start / $\mu$s")
    ax.set_ylabel("Shifted signal / mV")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")

    # Figure 76b: log-linear scale, exponential/lifetime fit.
    ax = axes[1]
    ax.semilogy(x_us[log_mask], y_mV[log_mask], ".", markersize=2, label="baseline-shifted data")
    if np.isfinite(tau2_s) and np.isfinite(slope2_per_s) and np.isfinite(intercept2_log):
        max_x = min(ANALYSIS_END_US, 200.0)
        fit_x_s = np.linspace(max(0.0, np.nanmin(t_rel_s[fit2_mask])), max_x * 1e-6, 250)
        fit_y_mV = np.exp(slope2_per_s * fit_x_s + intercept2_log) * 1e3
        ax.semilogy(fit_x_s * 1e6, fit_y_mV, linewidth=2, label=f"initial log-slope, tau = {tau2_s*1e6:.2f} $\mu$s")

    ax.set_xlabel("Time after decay start / $\mu$s")
    ax.set_ylabel("Shifted signal / mV")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend(loc="best")

    fig.suptitle(path.name)
    fig.savefig(figure_path, dpi=200)

    if not SHOW_PLOTS:
        plt.close(fig)


# ---------------------------------------------------------------------------
# MAIN PROGRAM
# ---------------------------------------------------------------------------

def find_input_files() -> list[Path]:
    if INPUT_FILES:
        return [Path(p) for p in INPUT_FILES]

    if DATA_FOLDER.exists():
        files: list[Path] = []
        for pattern in ("*.zip", "*.csv", "*.csv.gz"):
            files.extend(DATA_FOLDER.glob(pattern))
        return sorted(set(files))

    # No command line needed: if the configured folder is not found, show a file picker.
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        selected = filedialog.askopenfilenames(
            title="Select Tektronix waveform CSV/ZIP files",
            filetypes=[
                ("Waveform files", "*.zip *.csv *.csv.gz"),
                ("ZIP files", "*.zip"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        root.destroy()
        return [Path(p) for p in selected]
    except Exception:
        return []


def write_results_csv(results: list[LifetimeResult], output_folder: Path) -> Path:
    out = output_folder / "lifetime_results.csv"
    fieldnames = list(LifetimeResult.__dataclass_fields__.keys())
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)
    return out


def main() -> None:
    files = find_input_files()
    if not files:
        print("No input files found. Check DATA_FOLDER or select files in the dialog.")
        return

    missing = [p for p in files if not p.exists()]
    if missing:
        print("These files do not exist:")
        for p in missing:
            print(f"  {p}")
        return

    if OUTPUT_FOLDER is not None:
        output_folder = Path(OUTPUT_FOLDER)
    elif DATA_FOLDER.exists():
        output_folder = DATA_FOLDER / "lifetime_results"
    else:
        output_folder = files[0].parent / "lifetime_results"
    output_folder.mkdir(parents=True, exist_ok=True)

    print("Lifetime evaluation")
    print(f"Output folder: {output_folder}")
    print(f"Files: {len(files)}")

    results: list[LifetimeResult] = []
    for path in files:
        try:
            results.append(analyse_one_file(path, output_folder))
        except Exception as exc:
            print(f"\nERROR while analysing {path.name}: {exc}")

    if results:
        table_path = write_results_csv(results, output_folder)
        print(f"\nSaved results table: {table_path}")
        print("\nSummary:")
        for r in results:
            print(
                f"  {r.file}: "
                f"tau1 = {r.tau_method1_us:.3f} us, "
                f"tau2 = {r.tau_method2_us:.3f} us"
            )

    if SHOW_PLOTS and results:
        plt.show()


if __name__ == "__main__":
    try:
        main()
    finally:
        if PAUSE_AT_END:
            try:
                input("\nFinished. Press Enter to close...")
            except EOFError:
                pass
