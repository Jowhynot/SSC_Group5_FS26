from __future__ import annotations

import csv
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Use a non-interactive backend so the script works by double-clicking and just saves plots.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

# Your data folder. The script will use this folder if it exists.
# If it does not exist, a folder-selection window will open.
DATA_DIR = Path(r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5 ToF IV R\G5")

# Where output files are written. This folder is created inside DATA_DIR.
OUTPUT_FOLDER_NAME = "ToF_results"

# Signal channel in the Tektronix CSV. For your files, CH1 contains the ToF signal.
SIGNAL_CHANNEL = "CH1"

# Optional bias-monitor channel. The script mainly gets the applied voltage from the file name,
# e.g. 180V_TOF.csv -> 180 V. This channel is only used as a fallback.
BIAS_CHANNEL = "CH4"

# Contact distance L between electrodes/contacts.
# Set this if you know it and do not want the popup, e.g. CONTACT_DISTANCE_M = 1.0e-3 for 1 mm.
# Leave as None to be asked for L in mm when the program starts.
CONTACT_DISTANCE_M: Optional[float] = None

# Process files with "dark" in the name too? The user asked for every CSV file, so True.
PROCESS_DARK_FILES = True

# Smoothing used for derivative detection and fitting. 40 ns works well for the uploaded files.
SMOOTHING_TIME_S = 40e-9

# Ignore the first few ns after the laser edge in the log-log fit to avoid fitting the switching spike.
FIT_MIN_TIME_S = 50e-9

# Maximum time after laser edge included in fitting. None means automatic.
FIT_MAX_TIME_S: Optional[float] = None

# Keep this True to open the output folder automatically at the end on Windows/macOS/Linux.
OPEN_OUTPUT_FOLDER_WHEN_DONE = True

# Manual log-log fitting windows.
# Adjust these until the fits follow the physically meaningful regions.
PLATEAU_FIT_RANGE_US = (0.08, 1.0)
TAIL_FIT_RANGE_US = (2.2, 3)
TAIL_PLOT_RANGE_US = (0.01, 3)     # just for drawing a longer green line
PLATEAU_PLOT_RANGE_US = (0.01, 2.0)

LINEAR_LEVEL_V = 0.04

# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------

@dataclass
class FitResult:
    filename: str
    applied_voltage_V: Optional[float]
    laser_arrival_time_s: Optional[float]
    transient_time_s: Optional[float]
    mobility_m2_per_Vs: Optional[float]
    mobility_cm2_per_Vs: Optional[float]
    plateau_slope_loglog: Optional[float]
    tail_slope_loglog: Optional[float]
    fit_intersection_signal: Optional[float]
    polarity: Optional[int]
    status: str
    note: str


# -----------------------------------------------------------------------------
# Small UI helpers - lets the script be double-clickable without command-line args.
# -----------------------------------------------------------------------------

def _try_messagebox(title: str, message: str, error: bool = False) -> None:
    """Show a GUI message if possible; otherwise print to console."""
    try:
        import tkinter as tk
        from tkinter import messagebox

        root = tk.Tk()
        root.withdraw()
        if error:
            messagebox.showerror(title, message)
        else:
            messagebox.showinfo(title, message)
        root.destroy()
    except Exception:
        print(f"\n{title}\n{'=' * len(title)}\n{message}\n")


def choose_data_dir(default_dir: Path) -> Path:
    if default_dir.exists() and default_dir.is_dir():
        return default_dir

    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        selected = filedialog.askdirectory(
            title="Select folder containing your ToF CSV files",
            initialdir=str(Path.home()),
        )
        root.destroy()
        if selected:
            return Path(selected)
    except Exception:
        pass

    # Console fallback.
    typed = input("Enter the folder containing the ToF CSV files: ").strip().strip('"')
    if typed:
        return Path(typed)
    return default_dir


def get_contact_distance_m() -> Optional[float]:
    """Return L in meters, using CONTACT_DISTANCE_M or a small popup prompt."""
    if CONTACT_DISTANCE_M is not None and CONTACT_DISTANCE_M > 0:
        return float(CONTACT_DISTANCE_M)

    prompt = (
        "Enter the distance L between contacts in millimetres.\n\n"
        "Mobility is calculated as mu = L^2 / (t_tr * |U|).\n"
        "Cancel leaves mobility blank but still calculates transient times."
    )

    try:
        import tkinter as tk
        from tkinter import simpledialog

        root = tk.Tk()
        root.withdraw()
        value_mm = simpledialog.askfloat(
            "ToF contact distance",
            prompt,
            minvalue=0.000001,
        )
        root.destroy()
        if value_mm is None:
            return None
        return float(value_mm) * 1e-3
    except Exception:
        typed = input(prompt + "\nL in mm, or blank to skip mobility: ").strip()
        if not typed:
            return None
        try:
            return float(typed.replace(",", ".")) * 1e-3
        except ValueError:
            return None


def open_output_folder(path: Path) -> None:
    if not OPEN_OUTPUT_FOLDER_WHEN_DONE:
        return
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            import subprocess
            subprocess.Popen(["open", str(path)])
        else:
            import subprocess
            subprocess.Popen(["xdg-open", str(path)])
    except Exception:
        pass


# -----------------------------------------------------------------------------
# CSV loading
# -----------------------------------------------------------------------------


def find_header_line(csv_path: Path) -> Tuple[int, List[str]]:
    """Find the line containing TIME,CH1,... and return zero-based line index plus columns."""
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for line_no, line in enumerate(f):
            stripped = line.strip().lstrip("\ufeff")
            if stripped.upper().startswith("TIME,") or stripped.upper() == "TIME":
                cols = [c.strip().strip('"') for c in stripped.split(",")]
                return line_no, cols
    raise ValueError("Could not find a Tektronix data header line starting with TIME,")


def load_tektronix_csv(csv_path: Path) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Load TIME and channel arrays from a Tektronix oscilloscope CSV."""
    header_line, columns = find_header_line(csv_path)
    data = np.loadtxt(csv_path, delimiter=",", skiprows=header_line + 1)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    channels: Dict[str, np.ndarray] = {}
    for i, name in enumerate(columns):
        if i < data.shape[1] and name:
            channels[name.strip()] = np.asarray(data[:, i], dtype=float)

    if "TIME" not in channels:
        raise ValueError("No TIME column found.")
    return channels["TIME"], channels


def parse_voltage_from_filename(path: Path) -> Optional[float]:
    """Extract applied voltage from names such as 180V_TOF.csv or -200V_TOF_dark.csv."""
    match = re.search(r"(?<![A-Za-z0-9])([+-]?\d+(?:[.,]\d+)?)\s*V(?![A-Za-z0-9])", path.stem, re.IGNORECASE)
    if match:
        return float(match.group(1).replace(",", "."))
    return None


# -----------------------------------------------------------------------------
# Numerical helpers
# -----------------------------------------------------------------------------

def odd_at_least(n: int, minimum: int = 5) -> int:
    n = max(int(n), minimum)
    if n % 2 == 0:
        n += 1
    return n


def smooth_edge(y: np.ndarray, window: int) -> np.ndarray:
    """Moving-average smoothing with edge padding, avoiding fake edge dips."""
    if window <= 1:
        return y.astype(float, copy=True)
    window = min(window, len(y) - 1 if (len(y) - 1) % 2 else len(y) - 2)
    window = odd_at_least(window, 3)
    pad = window // 2
    padded = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(padded, kernel, mode="valid")


def robust_std(y: np.ndarray) -> float:
    if len(y) < 5:
        return float(np.std(y)) if len(y) else float("nan")
    med = np.median(y)
    mad = np.median(np.abs(y - med))
    if mad > 0:
        return float(1.4826 * mad)
    return float(np.std(y))


def downsample(x: np.ndarray, y: np.ndarray, max_points: int = 8000) -> Tuple[np.ndarray, np.ndarray]:
    if len(x) <= max_points:
        return x, y
    step = int(math.ceil(len(x) / max_points))
    return x[::step], y[::step]


def linear_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    coeff = np.polyfit(x, y, deg=1)
    return float(coeff[0]), float(coeff[1])

def falling_edge_level_time(
    tau_s: np.ndarray,
    signal: np.ndarray,
    level_v: float,
) -> Optional[float]:
    """
    Return the time after laser arrival where the falling signal crosses level_v.

    The function searches after the signal maximum and linearly interpolates
    between the two points around the crossing.
    """

    valid = (
        np.isfinite(tau_s)
        & np.isfinite(signal)
        & (tau_s > 0)
    )

    if np.count_nonzero(valid) < 5:
        return None

    tau = tau_s[valid]
    sig = signal[valid]

    order = np.argsort(tau)
    tau = tau[order]
    sig = sig[order]

    peak_idx = int(np.argmax(sig))

    for i in range(peak_idx + 1, len(sig)):
        y1 = sig[i - 1]
        y2 = sig[i]

        if y1 >= level_v and y2 <= level_v:
            x1 = tau[i - 1]
            x2 = tau[i]

            if y2 == y1:
                return float(x2)

            # Linear interpolation to get the exact crossing time
            frac = (level_v - y1) / (y2 - y1)
            return float(x1 + frac * (x2 - x1))

    return None

def fit_loglog_region(
    tau_s: np.ndarray,
    signal: np.ndarray,
    fit_range_us: Tuple[float, float],
    noise_level: float,
    min_points: int = 10,
) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """
    Fit log10(signal) = m * log10(time) + b
    using all valid data points inside a chosen time window.
    """

    lo_s = fit_range_us[0] * 1e-6
    hi_s = fit_range_us[1] * 1e-6

    mask = (
        np.isfinite(tau_s)
        & np.isfinite(signal)
        & (tau_s >= lo_s)
        & (tau_s <= hi_s)
        & (signal > max(5 * noise_level, 0))
    )

    if np.count_nonzero(mask) < min_points:
        raise ValueError(
            f"Too few valid points in fit range {fit_range_us[0]}-{fit_range_us[1]} us"
        )

    x = np.log10(tau_s[mask])
    y = np.log10(signal[mask])

    m, b = np.polyfit(x, y, deg=1)

    return float(m), float(b), tau_s[mask], signal[mask]


def fallback_half_height_time(tau_s: np.ndarray, signal: np.ndarray) -> Optional[float]:
    """Fallback if two-line log-log fit fails: first time after peak where signal drops to half of early level."""
    m = (tau_s > FIT_MIN_TIME_S) & np.isfinite(signal) & (signal > 0)
    if np.count_nonzero(m) < 20:
        return None
    tt = tau_s[m]
    yy = signal[m]
    max_index = int(np.argmax(yy))
    early_end = min(len(yy), max(max_index + 10, int(0.1 * len(yy))))
    early_level = float(np.percentile(yy[max_index:early_end], 75))
    if not np.isfinite(early_level) or early_level <= 0:
        return None
    half = early_level / 2.0
    after = np.where((np.arange(len(yy)) > max_index) & (yy <= half))[0]
    if len(after) == 0:
        return None
    return float(tt[int(after[0])])


# -----------------------------------------------------------------------------
# ToF analysis and plotting
# -----------------------------------------------------------------------------

def analyze_one_file(csv_path: Path, output_dir: Path, contact_distance_m: Optional[float]) -> FitResult:
    filename = csv_path.name
    applied_voltage = parse_voltage_from_filename(csv_path)

    try:
        t, columns = load_tektronix_csv(csv_path)

        if SIGNAL_CHANNEL not in columns:
            raise ValueError(f"Column {SIGNAL_CHANNEL!r} was not found. Available columns: {list(columns)}")

        v = np.asarray(columns[SIGNAL_CHANNEL], dtype=float)

        if len(t) < 100:
            raise ValueError("Not enough data points.")

        if applied_voltage is None and BIAS_CHANNEL in columns:
            bias = np.asarray(columns[BIAS_CHANNEL], dtype=float)
            applied_voltage = float(np.nanmedian(bias))

        dt = float(np.nanmedian(np.diff(t)))
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("TIME column is not strictly increasing.")

        smooth_points = odd_at_least(int(round(SMOOTHING_TIME_S / dt)), 5)
        smooth_points = min(smooth_points, odd_at_least(len(v) // 30, 5))
        v_smooth = smooth_edge(v, smooth_points)

        derivative = np.gradient(v_smooth, t)

        if float(np.nanmin(t)) < 0.0 < float(np.nanmax(t)):
            search_mask = (t > -2e-6) & (t < 2e-6)
            if np.count_nonzero(search_mask) < 50:
                search_mask = np.ones_like(t, dtype=bool)
        else:
            search_mask = np.ones_like(t, dtype=bool)

        margin = max(5, len(t) // 20)
        search_mask[:margin] = False
        search_mask[-margin:] = False

        search_indices = np.where(search_mask & np.isfinite(derivative))[0]
        if len(search_indices) == 0:
            raise ValueError("Could not search for the laser edge.")

        edge_idx = int(search_indices[np.argmax(np.abs(derivative[search_indices]))])
        laser_time = float(t[edge_idx])
        tau = t - laser_time

        pre_mask = t < (laser_time - 1.0e-6)
        late_start = laser_time + max(8e-6, 0.50 * (float(t[-1]) - laser_time))
        post_mask = t > late_start
        baseline_mask = pre_mask | post_mask

        if np.count_nonzero(baseline_mask) < 50:
            baseline_mask = pre_mask | (t > laser_time + 0.70 * (float(t[-1]) - laser_time))

        if np.count_nonzero(baseline_mask) >= 20:
            base_m, base_b = linear_fit(t[baseline_mask], v_smooth[baseline_mask])
            baseline = base_m * t + base_b
        else:
            baseline = np.full_like(
                v_smooth,
                float(np.median(v_smooth[: max(10, len(v_smooth) // 20)]))
            )

        corrected = v_smooth - baseline

        after_for_polarity = (
            (tau > FIT_MIN_TIME_S)
            & (tau < min(8e-6, 0.45 * (float(t[-1]) - laser_time)))
        )

        if np.count_nonzero(after_for_polarity) < 20:
            after_for_polarity = tau > FIT_MIN_TIME_S

        high = float(np.nanpercentile(corrected[after_for_polarity], 99))
        low = float(np.nanpercentile(corrected[after_for_polarity], 1))
        polarity = 1 if abs(high) >= abs(low) else -1
        signal = polarity * corrected

        noise_mask = pre_mask
        if np.count_nonzero(noise_mask) < 20:
            noise_mask = tau < -0.5e-6

        noise_level = robust_std(signal[noise_mask]) if np.count_nonzero(noise_mask) else robust_std(signal[: len(signal) // 10])
        max_signal = float(np.nanpercentile(signal[after_for_polarity], 99))

        note = ""

        try:
            m_plateau, b_plateau, plateau_tau, plateau_sig = fit_loglog_region(
                tau,
                signal,
                PLATEAU_FIT_RANGE_US,
                noise_level=noise_level,
            )

            m_tail, b_tail, tail_tau, tail_sig = fit_loglog_region(
                tau,
                signal,
                TAIL_FIT_RANGE_US,
                noise_level=noise_level,
            )

            denom = m_plateau - m_tail
            if abs(denom) < 1e-12:
                raise ValueError("Plateau and tail fits are nearly parallel.")

            log_t_tr = (b_tail - b_plateau) / denom
            t_tr = 10.0 ** log_t_tr
            y_intersection = 10.0 ** (m_plateau * log_t_tr + b_plateau)

            status = "ok"

        except Exception as fit_error:
            t_half = fallback_half_height_time(tau, signal)
            if t_half is None:
                raise fit_error

            t_tr = t_half
            m_plateau = b_plateau = m_tail = b_tail = y_intersection = float("nan")
            status = "fallback_half_height"
            note = f"Manual log-log fit failed; used half-height fallback. Reason: {fit_error}"

        mu_m2 = None
        mu_cm2 = None

        if (
            contact_distance_m is not None
            and applied_voltage is not None
            and abs(applied_voltage) > 0
            and t_tr is not None
            and np.isfinite(t_tr)
            and t_tr > 0
        ):
            mu_m2 = (contact_distance_m ** 2) / (float(t_tr) * abs(float(applied_voltage)))
            mu_cm2 = mu_m2 * 1.0e4

        make_figure73_plot(
            csv_path=csv_path,
            output_dir=output_dir,
            tau_s=tau,
            raw_v=v,
            signal=signal,
            transient_time_s=t_tr,
            plateau=(m_plateau, b_plateau),
            tail=(m_tail, b_tail),
            status=status,
        )

        return FitResult(
            filename=filename,
            applied_voltage_V=applied_voltage,
            laser_arrival_time_s=laser_time,
            transient_time_s=float(t_tr),
            mobility_m2_per_Vs=mu_m2,
            mobility_cm2_per_Vs=mu_cm2,
            plateau_slope_loglog=float(m_plateau) if np.isfinite(m_plateau) else None,
            tail_slope_loglog=float(m_tail) if np.isfinite(m_tail) else None,
            fit_intersection_signal=float(y_intersection) if np.isfinite(y_intersection) else None,
            polarity=polarity,
            status=status,
            note=note,
        )

    except Exception as exc:
        return FitResult(
            filename=filename,
            applied_voltage_V=applied_voltage,
            laser_arrival_time_s=None,
            transient_time_s=None,
            mobility_m2_per_Vs=None,
            mobility_cm2_per_Vs=None,
            plateau_slope_loglog=None,
            tail_slope_loglog=None,
            fit_intersection_signal=None,
            polarity=None,
            status="error",
            note=str(exc),
        )

def make_figure73_plot(
    csv_path: Path,
    output_dir: Path,
    tau_s: np.ndarray,
    raw_v: np.ndarray,
    signal: np.ndarray,
    transient_time_s: float,
    plateau: Tuple[float, float],
    tail: Tuple[float, float],
    status: str,
) -> None:
    figures_dir = output_dir / "figure_73_plots"
    figures_dir.mkdir(parents=True, exist_ok=True)
    safe_name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", csv_path.stem)
    out_png = figures_dir / f"{safe_name}.png"

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle(f"{csv_path.name}")

    # (a) Linear plot.
    x_us, y_sig = downsample(tau_s * 1e6, signal, max_points=10000)
    axes[0].plot(x_us, y_sig, linewidth=1.0)

    t_tr_A_s = falling_edge_level_time(
        tau_s=tau_s,
        signal=signal,
        level_v=LINEAR_LEVEL_V,
    )

    axes[0].axhline(
        LINEAR_LEVEL_V,
        linestyle=":",
        linewidth=1.2,
        label=f"{LINEAR_LEVEL_V:.2f} V level",
    )

    if t_tr_A_s is not None:
        axes[0].axvline(
            t_tr_A_s * 1e6,
            linestyle="--",
            linewidth=1.2,
            label=f"t_tr_A = {t_tr_A_s * 1e6:.2f} us",
        )
    else:
        axes[0].axvline(
            transient_time_s * 1e6,
            linestyle="--",
            linewidth=1.2,
            label="t_tr_A not found",
        )

    axes[0].axhline(0.0, linewidth=0.8)
    axes[0].set_xlabel("Time after laser arrival (μs)")
    axes[0].set_ylabel(f"Baseline-corrected {SIGNAL_CHANNEL} signal (V)")
    axes[0].grid(True, alpha=0.5)
    axes[0].legend(loc="best", fontsize=8)

    # (b) Log-log plot with line fits.
    fit_max = FIT_MAX_TIME_S
    if fit_max is None:
        fit_max = min(20e-6, 0.95 * float(np.nanmax(tau_s)))
    valid = (tau_s > FIT_MIN_TIME_S) & (tau_s < fit_max) & (signal > 0) & np.isfinite(signal)
    x_log_s = tau_s[valid]
    y_log = signal[valid]
    x_plot_us, y_plot = downsample(x_log_s * 1e6, y_log, max_points=12000)
    axes[1].loglog(x_plot_us, y_plot, linewidth=1.0, alpha=0.8, label="signal")

    m_plateau, b_plateau = plateau
    m_tail, b_tail = tail
    if np.isfinite(m_plateau) and np.isfinite(b_plateau) and np.isfinite(m_tail) and np.isfinite(b_tail):
        plateau_x_s = np.logspace(
            math.log10(PLATEAU_PLOT_RANGE_US[0] * 1e-6),
            math.log10(PLATEAU_PLOT_RANGE_US[1] * 1e-6),
            100,
        )
        
        tail_x_s = np.logspace(
            math.log10(TAIL_PLOT_RANGE_US[0] * 1e-6),
            math.log10(TAIL_PLOT_RANGE_US[1] * 1e-6),
            100,
        )
        
        axes[1].loglog(
            plateau_x_s * 1e6,
            10.0 ** (m_plateau * np.log10(plateau_x_s) + b_plateau),
            "--",
            linewidth=1.2,
            label="plateau fit",
        )
        
        axes[1].loglog(
            tail_x_s * 1e6,
            10.0 ** (m_tail * np.log10(tail_x_s) + b_tail),
            "--",
            linewidth=1.2,
            label="tail fit",
        )
    axes[1].axvline(
        transient_time_s * 1e6,
        linestyle="--",
        linewidth=1.2,
        label=f"t_tr = {transient_time_s * 1e6:.2f} µs",
    )
    axes[1].set_xlabel("Time after laser arrival / μs")
    axes[1].set_ylabel(f"Baseline-corrected {SIGNAL_CHANNEL} signal / V")
    axes[1].grid(True, which="both", alpha=0.5)
    axes[1].legend(loc="best", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def make_mobility_plot(results: List[FitResult], output_dir: Path) -> None:
    rows = [r for r in results if r.mobility_cm2_per_Vs is not None and r.applied_voltage_V is not None]
    if not rows:
        return
    rows.sort(key=lambda r: float(r.applied_voltage_V))
    volts = np.asarray([r.applied_voltage_V for r in rows], dtype=float)
    mob = np.asarray([r.mobility_cm2_per_Vs for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(volts, mob, marker="o", linestyle="-")
    for v, m, r in zip(volts, mob, rows):
        if "dark" in r.filename.lower():
            ax.annotate("dark", (v, m), textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax.set_xlabel("U / V ")
    ax.set_ylabel(r"Mobility $\mu$ / cm$^2$ V$^{-1}$ s$^{-1}$")
    ax.set_title("Charge-carrier mobility vs applied voltage")
    ax.grid(True, alpha=0.5)
    fig.tight_layout()
    fig.savefig(output_dir / "mobility_vs_voltage.png", dpi=200)
    plt.close(fig)


def write_summary_csv(results: List[FitResult], output_dir: Path, contact_distance_m: Optional[float]) -> Path:
    out_csv = output_dir / "tof_summary.csv"
    fields = [
        "filename",
        "applied_voltage_V",
        "laser_arrival_time_us",
        "transient_time_us",
        "contact_distance_m",
        "mobility_m2_per_Vs",
        "mobility_cm2_per_Vs",
        "plateau_slope_loglog",
        "tail_slope_loglog",
        "fit_intersection_signal_V",
        "polarity",
        "status",
        "note",
    ]

    def fmt(value: Optional[float]) -> str:
        if value is None:
            return ""
        try:
            if not np.isfinite(value):
                return ""
        except TypeError:
            pass
        return f"{value:.10g}"

    with out_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "filename": r.filename,
                "applied_voltage_V": fmt(r.applied_voltage_V),
                "laser_arrival_time_us": fmt(None if r.laser_arrival_time_s is None else r.laser_arrival_time_s * 1e6),
                "transient_time_us": fmt(None if r.transient_time_s is None else r.transient_time_s * 1e6),
                "contact_distance_m": fmt(contact_distance_m),
                "mobility_m2_per_Vs": fmt(r.mobility_m2_per_Vs),
                "mobility_cm2_per_Vs": fmt(r.mobility_cm2_per_Vs),
                "plateau_slope_loglog": fmt(r.plateau_slope_loglog),
                "tail_slope_loglog": fmt(r.tail_slope_loglog),
                "fit_intersection_signal_V": fmt(r.fit_intersection_signal),
                "polarity": "" if r.polarity is None else str(r.polarity),
                "status": r.status,
                "note": r.note,
            })
    return out_csv


# -----------------------------------------------------------------------------
# Main program
# -----------------------------------------------------------------------------

def main() -> None:
    data_dir = choose_data_dir(DATA_DIR)
    if not data_dir.exists() or not data_dir.is_dir():
        _try_messagebox("ToF analysis error", f"Data folder not found:\n{data_dir}", error=True)
        return

    csv_files = sorted(data_dir.glob("*.csv"))
    if not PROCESS_DARK_FILES:
        csv_files = [p for p in csv_files if "dark" not in p.name.lower()]

    if not csv_files:
        _try_messagebox("ToF analysis error", f"No CSV files found in:\n{data_dir}", error=True)
        return

    contact_distance_m = get_contact_distance_m()
    output_dir = data_dir / OUTPUT_FOLDER_NAME
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "tof_run_log.txt"
    results: List[FitResult] = []

    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"Data folder: {data_dir}\n")
        log.write(f"Output folder: {output_dir}\n")
        log.write(f"Contact distance L: {contact_distance_m} m\n")
        log.write(f"Number of CSV files: {len(csv_files)}\n\n")

        for i, csv_path in enumerate(csv_files, start=1):
            print(f"[{i}/{len(csv_files)}] Processing {csv_path.name} ...")
            log.write(f"[{i}/{len(csv_files)}] {csv_path.name}\n")
            try:
                result = analyze_one_file(csv_path, output_dir, contact_distance_m)
                results.append(result)
                log.write(f"  status: {result.status}\n")
                if result.note:
                    log.write(f"  note: {result.note}\n")
            except Exception:
                tb = traceback.format_exc()
                log.write(tb + "\n")
                results.append(FitResult(
                    filename=csv_path.name,
                    applied_voltage_V=parse_voltage_from_filename(csv_path),
                    laser_arrival_time_s=None,
                    transient_time_s=None,
                    mobility_m2_per_Vs=None,
                    mobility_cm2_per_Vs=None,
                    plateau_slope_loglog=None,
                    tail_slope_loglog=None,
                    fit_intersection_signal=None,
                    polarity=None,
                    status="error",
                    note=tb.splitlines()[-1] if tb else "unknown error",
                ))

    summary_csv = write_summary_csv(results, output_dir, contact_distance_m)
    make_mobility_plot(results, output_dir)

    ok_count = sum(1 for r in results if r.status in {"ok", "fallback_half_height"})
    error_count = len(results) - ok_count
    mobility_note = ""
    if contact_distance_m is None:
        mobility_note = "\n\nMobility was left blank because no contact distance L was entered."

    message = (
        f"ToF analysis finished.\n\n"
        f"Processed files: {len(results)}\n"
        f"Successful analyses: {ok_count}\n"
        f"Errors: {error_count}\n\n"
        f"Summary table:\n{summary_csv}\n\n"
        f"Figure-73 plots:\n{output_dir / 'figure_73_plots'}\n\n"
        f"Run log:\n{log_path}"
        f"{mobility_note}"
    )
    _try_messagebox("ToF analysis complete", message, error=(error_count > 0 and ok_count == 0))
    open_output_folder(output_dir)


if __name__ == "__main__":
    main()
