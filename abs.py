from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import t as student_t
except Exception:  # pragma: no cover
    student_t = None


# ----------------------------- user settings -----------------------------

PLOT_ENERGY_MIN_EV = 2.1
PLOT_ENERGY_MAX_EV = 2.4

DATA_DIR = Path(r"C:\Users\jsjan\Downloads\SSC Data&Codes\Group 5_abs_photolumin\Group 5")

# Direct allowed band gap: r = 1/2 -> plot (alpha*hnu)^2 against hnu.
TRANSITION_R = 0.5

CONFIDENCE_LEVEL = 0.95
MIN_FIT_POINTS = 12

# Used only when thickness cannot be inferred from filename/sample name.
# The Tauc intercept and Urbach slope are not changed by multiplying alpha by a
# constant, but alpha's physical cm^-1 scale does need the sample thickness.
DEFAULT_THICKNESS_CM = 1.0

# Common thickness rules for your files.
# 5 µm = 5e-4 cm; 2 mm = 0.2 cm.
THICKNESS_RULES: List[Tuple[str, float]] = [
    # Put single-crystal rules first because some exported CSV filenames can be
    # misleading copies, while the sample name inside the CSV is usually reliable.
    (r"single\s*crystal|singlecrystal|2\s*mm|2mmsw|2mm", 0.2),
    (r"5\s*(ym|um|µm|u?m)|5ym|5um", 5e-4),
]

# Optional manual fit windows in eV.
# Keys may be the exact filename, file stem, or sample name found in the CSV.
# Example:
# MANUAL_URBACH_RANGES = {
#     "G5_5ym_film": (2.30, 2.36),
#     "singlecrystal_absorption_sw_2nm_300_650nm_(1).csv": (2.22, 2.28),
# }
MANUAL_URBACH_RANGES: Dict[str, Tuple[float, float]] = {}
MANUAL_TAUC_RANGES: Dict[str, Tuple[float, float]] = {}

# Automatic fit-window thresholds. Edit only if the automatic windows are poor.
URBACH_NORM_ABSORBANCE_RANGE = (0.02, 0.45)   # exponential tail / onset
TAUC_NORM_ABSORBANCE_RANGE = (0.12, 0.88)     # rising absorption edge

HC_EV_NM = 1239.841984  # h*c in eV nm


# ----------------------------- data classes ------------------------------

@dataclass
class LinearFit:
    slope: float
    intercept: float
    slope_ci: Tuple[float, float]
    intercept_ci: Tuple[float, float]
    cov: np.ndarray
    r2: float
    n: int
    dof: int
    x_min: float
    x_max: float


@dataclass
class SampleResult:
    filename: str
    sample: str
    safe_name: str
    thickness_cm: float
    wavelength_nm: np.ndarray
    energy_eV: np.ndarray
    absorbance: np.ndarray
    alpha_cm: np.ndarray
    urbach_fit: LinearFit
    tauc_fit: LinearFit
    urbach_energy_eV: float
    urbach_energy_ci_eV: Tuple[float, float]
    tauc_gap_eV: float
    tauc_gap_ci_eV: Tuple[float, float]
    transition_r: float
    warnings: List[str]


# ----------------------------- utilities ---------------------------------

def sanitize_filename(name: str, max_len: int = 120) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return safe[:max_len] if safe else "sample"


def try_tcrit(confidence: float, dof: int) -> float:
    if student_t is not None and dof > 0:
        return float(student_t.ppf(0.5 + confidence / 2.0, dof))
    return 1.96


def parse_sample_name(path: Path) -> str:
    """Try to read Shimadzu's sample name from the first line; fallback to stem."""
    try:
        first = path.read_text(errors="replace").splitlines()[0]
        parts = [p.strip().strip('"') for p in first.split(",")]
        for p in parts:
            if p:
                return p
    except Exception:
        pass
    return path.stem


def read_absorption_csv(path: Path) -> Tuple[str, pd.DataFrame]:
    """Read CSV and return sample name plus numeric wavelength/absorbance table."""
    sample = parse_sample_name(path)

    # sep=None lets pandas sniff comma/semicolon/tab delimiters.
    try:
        raw = pd.read_csv(path, header=None, sep=None, engine="python")
    except Exception:
        raw = pd.read_csv(path, header=None)

    if raw.shape[1] < 2:
        raise ValueError("CSV must contain at least two columns: wavelength and absorbance.")

    wavelength = pd.to_numeric(raw.iloc[:, 0], errors="coerce")
    absorbance = pd.to_numeric(raw.iloc[:, 1], errors="coerce")

    mask = wavelength.notna() & absorbance.notna()
    mask &= wavelength.between(100, 2500)  # nm sanity range
    df = pd.DataFrame({
        "wavelength_nm": wavelength[mask].astype(float),
        "absorbance": absorbance[mask].astype(float),
    })

    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    df = df[df["absorbance"] > 0]
    df = df.drop_duplicates(subset=["wavelength_nm"]).sort_values("wavelength_nm")

    if len(df) < 20:
        raise ValueError(f"Only {len(df)} numeric data points found; check CSV format.")

    return sample, df.reset_index(drop=True)


def infer_thickness_cm(path: Path, sample: str) -> Tuple[float, str]:
    # Prefer the sample name stored inside the CSV, then fall back to filename.
    # This matters when a copied/exported file has a misleading filename.
    keys = [
        sample.lower().replace("µ", "u"),
        f"{path.name} {path.stem}".lower().replace("µ", "u"),
    ]
    for key in keys:
        for pattern, thickness in THICKNESS_RULES:
            if re.search(pattern, key, flags=re.IGNORECASE):
                return thickness, ""
    return DEFAULT_THICKNESS_CM, (
        f"Thickness not inferred; used DEFAULT_THICKNESS_CM={DEFAULT_THICKNESS_CM:g} cm. "
        "Edit THICKNESS_RULES or DEFAULT_THICKNESS_CM if needed."
    )


def get_manual_range(
    mapping: Dict[str, Tuple[float, float]], path: Path, sample: str
) -> Optional[Tuple[float, float]]:
    for key in (path.name, path.stem, sample):
        if key in mapping:
            lo, hi = mapping[key]
            return (min(lo, hi), max(lo, hi))
    return None


def linear_fit_ci(x: np.ndarray, y: np.ndarray, confidence: float) -> LinearFit:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 3:
        raise ValueError("Need at least 3 points for a confidence interval.")

    X = np.column_stack([x, np.ones_like(x)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    slope, intercept = float(beta[0]), float(beta[1])

    y_hat = X @ beta
    residuals = y - y_hat
    n = len(x)
    dof = n - 2
    sse = float(np.sum(residuals ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else float("nan")

    mse = sse / dof if dof > 0 else float("nan")
    xtx_inv = np.linalg.pinv(X.T @ X)
    cov = mse * xtx_inv

    tcrit = try_tcrit(confidence, dof)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    slope_ci = (slope - tcrit * se[0], slope + tcrit * se[0])
    intercept_ci = (intercept - tcrit * se[1], intercept + tcrit * se[1])

    return LinearFit(
        slope=slope,
        intercept=intercept,
        slope_ci=slope_ci,
        intercept_ci=intercept_ci,
        cov=cov,
        r2=r2,
        n=n,
        dof=dof,
        x_min=float(np.min(x)),
        x_max=float(np.max(x)),
    )


def transformed_ci(
    value: float,
    gradient: np.ndarray,
    cov: np.ndarray,
    dof: int,
    confidence: float,
) -> Tuple[float, float]:
    variance = float(gradient @ cov @ gradient.T)
    se = math.sqrt(max(variance, 0.0))
    tcrit = try_tcrit(confidence, dof)
    return (value - tcrit * se, value + tcrit * se)


def find_segments(mask: np.ndarray, min_points: int, max_gap: int = 3) -> List[Tuple[int, int]]:
    inds = np.flatnonzero(mask)
    if len(inds) == 0:
        return []

    segments: List[Tuple[int, int]] = []
    start = int(inds[0])
    prev = int(inds[0])

    for idx in inds[1:]:
        idx = int(idx)
        if idx > prev + max_gap:
            if prev - start + 1 >= min_points:
                segments.append((start, prev))
            start = idx
        prev = idx

    if prev - start + 1 >= min_points:
        segments.append((start, prev))

    return segments


def score_linear_window(fit: LinearFit, y_window: np.ndarray) -> float:
    y_range = float(np.nanmax(y_window) - np.nanmin(y_window))
    return fit.r2 + 0.03 * math.log(max(fit.n, 2)) + 0.05 * min(y_range, 2.0) / 2.0


def choose_urbach_range(
    energy_eV: np.ndarray,
    alpha_cm: np.ndarray,
    absorbance: np.ndarray,
    manual_range: Optional[Tuple[float, float]],
) -> Tuple[np.ndarray, str]:
    x = energy_eV
    y = np.log(alpha_cm)

    if manual_range is not None:
        lo, hi = manual_range
        mask = (x >= lo) & (x <= hi)
        if np.count_nonzero(mask) < MIN_FIT_POINTS:
            raise ValueError(f"Manual Urbach range {manual_range} contains too few points.")
        return mask, "manual"

    lo_a, hi_a = np.nanpercentile(absorbance, [5, 85])
    if hi_a <= lo_a:
        norm = absorbance / np.nanmax(absorbance)
    else:
        norm = (absorbance - lo_a) / (hi_a - lo_a)

    lo_n, hi_n = URBACH_NORM_ABSORBANCE_RANGE
    allowed = (norm >= lo_n) & (norm <= hi_n)

    best_mask = None
    best_score = -np.inf
    widths_eV = [0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25]

    for width in widths_eV:
        for i in range(len(x)):
            j = int(np.searchsorted(x, x[i] + width, side="right"))
            if j - i < MIN_FIT_POINTS:
                continue
            if np.mean(allowed[i:j]) < 0.80:
                continue

            xx = x[i:j]
            yy = y[i:j]
            if np.nanmax(yy) - np.nanmin(yy) < 0.10:
                continue

            fit = linear_fit_ci(xx, yy, CONFIDENCE_LEVEL)
            if fit.slope <= 0:
                continue

            eu = 1.0 / fit.slope
            if not (0.001 <= eu <= 1.0):  # 1 meV to 1000 meV sanity range
                continue

            score = score_linear_window(fit, yy)
            if score > best_score:
                best_score = score
                mask = np.zeros_like(x, dtype=bool)
                mask[i:j] = True
                best_mask = mask

    if best_mask is None:
        # Conservative fallback: use the first segment in the tail/onset region.
        segments = find_segments(allowed, MIN_FIT_POINTS)
        if not segments:
            raise ValueError("Could not find an automatic Urbach fit window.")
        s, e = segments[0]
        best_mask = np.zeros_like(x, dtype=bool)
        best_mask[s:e + 1] = True
        return best_mask, "auto_fallback"

    return best_mask, "auto"


def choose_tauc_range(
    energy_eV: np.ndarray,
    alpha_cm: np.ndarray,
    absorbance: np.ndarray,
    transition_r: float,
    manual_range: Optional[Tuple[float, float]],
) -> Tuple[np.ndarray, str]:
    x = energy_eV
    y = (alpha_cm * x) ** (1.0 / transition_r)

    if manual_range is not None:
        lo, hi = manual_range
        mask = (x >= lo) & (x <= hi)
        if np.count_nonzero(mask) < MIN_FIT_POINTS:
            raise ValueError(f"Manual Tauc range {manual_range} contains too few points.")
        return mask, "manual"

    lo_a, hi_a = np.nanpercentile(absorbance, [5, 85])
    if hi_a <= lo_a:
        norm = absorbance / np.nanmax(absorbance)
    else:
        norm = (absorbance - lo_a) / (hi_a - lo_a)

    lo_n, hi_n = TAUC_NORM_ABSORBANCE_RANGE
    edge_allowed = (norm >= lo_n) & (norm <= hi_n)

    segments = find_segments(edge_allowed, MIN_FIT_POINTS)
    if segments:
        # Use the first rising-edge segment in ascending energy. This avoids
        # selecting high-energy linear-looking plateaus as the Tauc region.
        s0, e0 = segments[0]
        search_start, search_end = s0, e0 + 1
    else:
        search_start, search_end = 0, len(x)

    # Normalize only for scoring; final fit uses the physical y values.
    y_search = y[search_start:search_end]
    scale = np.nanmax(np.abs(y_search))
    if not np.isfinite(scale) or scale == 0:
        scale = 1.0
    y_norm = y / scale

    best_mask = None
    best_score = -np.inf
    widths_eV = [0.025, 0.035, 0.05, 0.07, 0.09, 0.12, 0.15]

    for width in widths_eV:
        for i in range(search_start, search_end):
            j = int(np.searchsorted(x, x[i] + width, side="right"))
            j = min(j, search_end)
            if j - i < MIN_FIT_POINTS:
                continue

            xx = x[i:j]
            yy = y_norm[i:j]
            if np.nanmax(yy) - np.nanmin(yy) < 0.03:
                continue

            fit = linear_fit_ci(xx, yy, CONFIDENCE_LEVEL)
            if fit.slope <= 0:
                continue

            eg = -fit.intercept / fit.slope
            # Tauc intercept should lie slightly below the chosen fit window.
            if not (x[0] - 0.2 <= eg <= xx[0] + 0.08):
                continue

            score = score_linear_window(fit, yy)
            if score > best_score:
                best_score = score
                mask = np.zeros_like(x, dtype=bool)
                mask[i:j] = True
                best_mask = mask

    if best_mask is None:
        if not segments:
            raise ValueError("Could not find an automatic Tauc fit window.")
        s, e = segments[0]
        best_mask = np.zeros_like(x, dtype=bool)
        best_mask[s:e + 1] = True
        return best_mask, "auto_fallback"

    return best_mask, "auto"


def process_file(path: Path, transition_r: float, confidence: float, out_dir: Path) -> SampleResult:
    sample, df = read_absorption_csv(path)
    thickness_cm, thickness_warning = infer_thickness_cm(path, sample)

    wavelength = df["wavelength_nm"].to_numpy(dtype=float)
    absorbance = df["absorbance"].to_numpy(dtype=float)
    energy = HC_EV_NM / wavelength
    alpha = 2.303 * absorbance / thickness_cm

    # Sort by increasing photon energy.
    order = np.argsort(energy)
    wavelength = wavelength[order]
    energy = energy[order]
    absorbance = absorbance[order]
    alpha = alpha[order]

    finite = (
        np.isfinite(wavelength)
        & np.isfinite(energy)
        & np.isfinite(absorbance)
        & np.isfinite(alpha)
        & (alpha > 0)
        & (absorbance > 0)
    )
    wavelength = wavelength[finite]
    energy = energy[finite]
    absorbance = absorbance[finite]
    alpha = alpha[finite]

    warnings: List[str] = []
    if thickness_warning:
        warnings.append(thickness_warning)

    if np.nanmax(absorbance) > 3.0:
        warnings.append(
            "Absorbance exceeds 3 for part of the spectrum; detector saturation/scattering may affect the fit."
        )

    urbach_manual = get_manual_range(MANUAL_URBACH_RANGES, path, sample)
    tauc_manual = get_manual_range(MANUAL_TAUC_RANGES, path, sample)

    urbach_mask, urbach_mode = choose_urbach_range(energy, alpha, absorbance, urbach_manual)
    tauc_mask, tauc_mode = choose_tauc_range(energy, alpha, absorbance, transition_r, tauc_manual)

    if urbach_mode != "manual":
        warnings.append(f"Urbach fit window selected automatically: {urbach_mode}. Inspect plot.")
    if tauc_mode != "manual":
        warnings.append(f"Tauc fit window selected automatically: {tauc_mode}. Inspect plot.")

    # Urbach: ln(alpha) = intercept + energy / E_U.
    urbach_fit = linear_fit_ci(energy[urbach_mask], np.log(alpha[urbach_mask]), confidence)
    urbach_energy = 1.0 / urbach_fit.slope
    urbach_ci = transformed_ci(
        urbach_energy,
        np.array([-1.0 / (urbach_fit.slope ** 2), 0.0]),
        urbach_fit.cov,
        urbach_fit.dof,
        confidence,
    )

    # Tauc: y = mE + b, Eg = -b/m.
    tauc_y = (alpha * energy) ** (1.0 / transition_r)
    tauc_fit = linear_fit_ci(energy[tauc_mask], tauc_y[tauc_mask], confidence)
    tauc_gap = -tauc_fit.intercept / tauc_fit.slope
    tauc_ci = transformed_ci(
        tauc_gap,
        np.array([tauc_fit.intercept / (tauc_fit.slope ** 2), -1.0 / tauc_fit.slope]),
        tauc_fit.cov,
        tauc_fit.dof,
        confidence,
    )

    safe_name = sanitize_filename(f"{path.stem}_{sample}")

    result = SampleResult(
        filename=path.name,
        sample=sample,
        safe_name=safe_name,
        thickness_cm=thickness_cm,
        wavelength_nm=wavelength,
        energy_eV=energy,
        absorbance=absorbance,
        alpha_cm=alpha,
        urbach_fit=urbach_fit,
        tauc_fit=tauc_fit,
        urbach_energy_eV=urbach_energy,
        urbach_energy_ci_eV=urbach_ci,
        tauc_gap_eV=tauc_gap,
        tauc_gap_ci_eV=tauc_ci,
        transition_r=transition_r,
        warnings=warnings,
    )

    plot_urbach(result, out_dir / f"{safe_name}_urbach.png")
    plot_tauc(result, out_dir / f"{safe_name}_tauc.png")

    return result


# ----------------------------- plotting ----------------------------------

def plot_urbach(result: SampleResult, out_path: Path) -> None:
    x = result.energy_eV
    alpha = result.alpha_cm
    fit = result.urbach_fit

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    plot_mask = (x >= PLOT_ENERGY_MIN_EV) & (x <= PLOT_ENERGY_MAX_EV)

    ax.semilogy(
        x[plot_mask],
        alpha[plot_mask],
        marker=".",
        linestyle="None",
        markersize=2.5,
        label=result.sample,
    )

    fit_line_min = max(fit.x_min, PLOT_ENERGY_MIN_EV)
    fit_line_max = min(fit.x_max, PLOT_ENERGY_MAX_EV)

    if fit_line_min < fit_line_max:
        fit_x = np.linspace(fit_line_min, fit_line_max, 200)
        fit_alpha = np.exp(fit.slope * fit_x + fit.intercept)
        ax.semilogy(fit_x, fit_alpha, linewidth=2.0, label="Urbach fit")

        ax.axvspan(fit_line_min, fit_line_max, alpha=0.12)

    eu_mev = 1000.0 * result.urbach_energy_eV
    ci_mev = (1000.0 * result.urbach_energy_ci_eV[0], 1000.0 * result.urbach_energy_ci_eV[1])
    text = (
        f"E$_U$ = {eu_mev:.1f} meV\n"
        f"95% CI [{ci_mev[0]:.1f}, {ci_mev[1]:.1f}] meV\n"
        f"R$^2$ = {fit.r2:.4f}"
    )
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top")

    ax.set_xlabel("Energy hν / eV")
    ax.set_ylabel(r"Absorption coefficient $\alpha$ / cm$^{-1}$")
    ax.set_xlim(PLOT_ENERGY_MIN_EV, PLOT_ENERGY_MAX_EV)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_tauc(result: SampleResult, out_path: Path) -> None:
    x = result.energy_eV
    y = (result.alpha_cm * x) ** (1.0 / result.transition_r)
    fit = result.tauc_fit

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    plot_mask = (x >= PLOT_ENERGY_MIN_EV) & (x <= PLOT_ENERGY_MAX_EV)

    ax.plot(
        x[plot_mask],
        y[plot_mask],
        marker=".",
        linestyle="None",
        markersize=2.5,
        label=result.sample,
    )

    # Plot line over fit window plus extrapolation to the intercept,
    # but only inside the displayed energy range.
    x_line_min = min(result.tauc_gap_eV, fit.x_min)
    x_line_max = fit.x_max

    fit_line_min = max(x_line_min, PLOT_ENERGY_MIN_EV)
    fit_line_max = min(x_line_max, PLOT_ENERGY_MAX_EV)

    if fit_line_min < fit_line_max:
        fit_x = np.linspace(fit_line_min, fit_line_max, 200)
        fit_y = fit.slope * fit_x + fit.intercept
        fit_line, = ax.plot(fit_x, fit_y, linewidth=2.0, label="Tauc fit")
    else:
        fit_line = None

    ax.axvspan(fit.x_min, fit.x_max, alpha=0.12)
    ax.axhline(0, linewidth=0.8)
    # Mark the intercept between the Tauc fit line and the baseline y=0.
    if PLOT_ENERGY_MIN_EV <= result.tauc_gap_eV <= PLOT_ENERGY_MAX_EV:
        marker_color = fit_line.get_color() if fit_line is not None else None

        ax.plot(
            result.tauc_gap_eV,
            0,
            marker="o",
            markersize=8,
            markerfacecolor="none",
            markeredgewidth=1.8,
            color=marker_color,
            linestyle="None",
            label="Tauc intercept",
        )

    eg = result.tauc_gap_eV
    ci = result.tauc_gap_ci_eV
    text = (
        f"E$_g$ = {eg:.4f} eV\n"
        f"95% CI [{ci[0]:.4f}, {ci[1]:.4f}] eV\n"
        f"R$^2$ = {fit.r2:.4f}"
    )
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top")

    exponent = 1.0 / result.transition_r
    if abs(exponent - round(exponent)) < 1e-10:
        exponent_label = f"{int(round(exponent))}"
    else:
        exponent_label = f"{exponent:.3g}"

    ax.set_xlabel("Energy hν / eV")
    ax.set_ylabel(rf"$(\alpha h\nu)^{{{exponent_label}}}$ / a.u.")
    ax.set_xlim(PLOT_ENERGY_MIN_EV, PLOT_ENERGY_MAX_EV)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)




# ----------------------------- output ------------------------------------

def summary_row(res: SampleResult) -> Dict[str, object]:
    return {
        "filename": res.filename,
        "sample": res.sample,
        "thickness_cm": res.thickness_cm,
        "transition_r": res.transition_r,
        "urbach_fit_E_min_eV": res.urbach_fit.x_min,
        "urbach_fit_E_max_eV": res.urbach_fit.x_max,
        "urbach_R2": res.urbach_fit.r2,
        "urbach_EU_eV": res.urbach_energy_eV,
        "urbach_EU_meV": 1000.0 * res.urbach_energy_eV,
        "urbach_EU_CI_low_meV": 1000.0 * res.urbach_energy_ci_eV[0],
        "urbach_EU_CI_high_meV": 1000.0 * res.urbach_energy_ci_eV[1],
        "tauc_fit_E_min_eV": res.tauc_fit.x_min,
        "tauc_fit_E_max_eV": res.tauc_fit.x_max,
        "tauc_R2": res.tauc_fit.r2,
        "tauc_Eg_eV": res.tauc_gap_eV,
        "tauc_Eg_CI_low_eV": res.tauc_gap_ci_eV[0],
        "tauc_Eg_CI_high_eV": res.tauc_gap_ci_eV[1],
        "warnings": " | ".join(res.warnings),
    }


def print_summary(summary: pd.DataFrame, out_dir: Path) -> None:
    cols = [
        "filename",
        "sample",
        "urbach_EU_meV",
        "urbach_EU_CI_low_meV",
        "urbach_EU_CI_high_meV",
        "tauc_Eg_eV",
        "tauc_Eg_CI_low_eV",
        "tauc_Eg_CI_high_eV",
        "urbach_R2",
        "tauc_R2",
    ]
    view = summary[cols].copy()

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 180)
    print("\nResults")
    print(view.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print(f"\nSaved figures and table to:\n  {out_dir}")


# ----------------------------- main --------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate Urbach energy and Tauc band gap from absorption CSV files.")
    parser.add_argument("--data-dir", type=str, default=str(DATA_DIR), help="Folder containing absorption CSV files.")
    parser.add_argument("--out-dir", type=str, default=None, help="Output folder. Default: DATA_DIR/abs_results")
    parser.add_argument("--transition-r", type=float, default=TRANSITION_R, help="Tauc transition exponent r. Direct allowed: 0.5.")
    parser.add_argument("--confidence", type=float, default=CONFIDENCE_LEVEL, help="Confidence level, e.g. 0.95.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser() if args.out_dir else data_dir / "abs_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(data_dir.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in: {data_dir}")

    results: List[SampleResult] = []
    errors: List[Tuple[str, str]] = []

    for path in csv_files:
        try:
            res = process_file(path, args.transition_r, args.confidence, out_dir)
            results.append(res)
        except Exception as exc:
            errors.append((path.name, str(exc)))
            print(f"WARNING: skipped {path.name}: {exc}")

    if not results:
        raise RuntimeError("No files were processed successfully.")

    summary = pd.DataFrame([summary_row(res) for res in results])
    summary.to_csv(out_dir / "abs_summary.csv", index=False)

    if errors:
        err = pd.DataFrame(errors, columns=["filename", "error"])
        err.to_csv(out_dir / "abs_errors.csv", index=False)

    print_summary(summary, out_dir)



if __name__ == "__main__":
    main()