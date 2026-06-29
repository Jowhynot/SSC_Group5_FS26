import csv
import math
import os
import re
import sys
import traceback
import zipfile
from pathlib import Path

# =============================================================================
# USER SETTINGS - EDIT THIS SECTION
# =============================================================================

# Your data folder. This can also be a .zip file path.
DATA_PATH = r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5 ToF IV R\G5"

# If DATA_PATH does not exist, the script will also try to find IV.zip or .dat
# files in the same folder as this script, then it will ask you to select a path.
RECURSIVE_SEARCH = True
OUTPUT_FOLDER_NAME = "IVR_results"

# Geometry must be in cm units to get rho in Ohm*cm:
#     A = contact area in cm^2
#     L = distance between contacts / crystal thickness in cm
#
# Fill the values after you measure the crystals.
# You may enter values for each contact/device key or only DEFAULT.
# Keys are extracted from filenames, e.g.:
#   IV_positive_gold_gold_dark__1.dat       -> gold_gold
#   IV_negative_gold_Bi_dark__1.dat         -> gold_Bi
#   IV_positive_Bi_gold_dark__1.dat         -> Bi_gold
#   IV_positive_FTO_right_way_dark__1.dat   -> FTO_right_way
#
# Example:
#   "gold_gold": {"area_cm2": 0.0314, "distance_cm": 0.10},
GEOMETRY_BY_CONTACT = {
    "DEFAULT":      {"area_cm2": 0.192, "distance_cm": 0.2},
    "gold_gold":    {"area_cm2": 0.192, "distance_cm": 0.2},
    "gold_Bi":      {"area_cm2": 0.192, "distance_cm": 0.2},
    "Bi_gold":      {"area_cm2": 0.192, "distance_cm": 0.2},
    "FTO":          {"area_cm2": 0.0064, "distance_cm": 0.01},
    "FTO_right_way": {"area_cm2": 0.0064, "distance_cm": 0.01},
}

# If you already know a manual linear/ohmic region, enter it here as |V| ranges.
# The script will then fit only points with vmin <= abs(V) <= vmax.
# If a key is not listed or is None, the script selects a linear region automatically.
# Example:
#   FIT_ABS_V_RANGE_BY_CONTACT = {"gold_gold": (0.0, 5.0), "FTO": (0.0, 0.5)}
FIT_ABS_V_RANGE_BY_CONTACT = {
    # "DEFAULT": (0.0, 5.0),
    # "gold_gold": (0.0, 5.0),
}

# Automatic linear-region selection settings. The script searches each sweep branch
# for the longest voltage window with R^2 >= AUTO_MIN_R2. If none is found, it uses
# the best available window and marks it in the warning column/report.
AUTO_MIN_POINTS = 6
AUTO_MIN_R2 = 0.990
AUTO_MIN_SPAN_FRACTION = 0.10

# Thresholds for warnings in positive/negative and hysteresis comparisons.
DIRECTION_RATIO_WARN = 2.0       # warn if max(R_pos, R_neg)/min(R_pos, R_neg) >= this
HYSTERESIS_RATIO_WARN = 2.0      # warn if two sweep branches differ by this factor
POOR_FIT_R2_WARN = 0.990

# Literature ranges for CsPbBr3 resistivity in Ohm*cm.
# Edit/add references to match the exact material/contact system used in your lab.
# The broad window below combines values from CsPbBr3 single-crystal reports:
# ~0.10 +/- 0.04 GOhm*cm, ~7e9 Ohm*cm, and up to ~2.5e10 Ohm*cm.
LITERATURE_RANGES_OHM_CM = [
    {"label": "CsPbBr3 solution-grown SC, dark, carbon contacts", "min": 6.0e7, "max": 1.4e8},
    {"label": "CsPbBr3 Bridgman SC, Ni/CsPbBr3/Ni ohmic contacts", "min": 7.0e9, "max": 7.0e9},
    {"label": "Broad CsPbBr3 single-crystal literature window", "min": 6.0e7, "max": 2.5e10},
]

MAKE_PLOTS = True
PAUSE_AT_END = os.environ.get("IVR_PAUSE", "1") != "0"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def circle_area_cm2(diameter_mm):
    """Optional helper: contact area for a circular contact from diameter in mm."""
    radius_cm = (diameter_mm / 10.0) / 2.0
    return math.pi * radius_cm * radius_cm


def rectangle_area_cm2(width_mm, height_mm):
    """Optional helper: contact area for a rectangular contact from mm dimensions."""
    return (width_mm / 10.0) * (height_mm / 10.0)


def safe_float(text):
    if text is None:
        return None
    text = str(text).strip().replace("D", "E").replace("d", "E")
    # Decimal comma support, but do not damage scientific notation with commas as separators.
    if text.count(",") == 1 and "." not in text:
        text = text.replace(",", ".")
    try:
        value = float(text)
        if math.isfinite(value):
            return value
    except ValueError:
        pass
    return None


def median(values):
    values = sorted([v for v in values if v is not None and math.isfinite(v)])
    if not values:
        return None
    n = len(values)
    mid = n // 2
    if n % 2:
        return values[mid]
    return 0.5 * (values[mid - 1] + values[mid])


def ratio_max_min(a, b):
    if a is None or b is None or a <= 0 or b <= 0:
        return None
    return max(a, b) / min(a, b)


def percent_difference(pos_value, neg_value):
    """Symmetric percent difference: 200*(pos-neg)/(pos+neg)."""
    if pos_value is None or neg_value is None:
        return None
    denom = pos_value + neg_value
    if denom == 0:
        return None
    return 200.0 * (pos_value - neg_value) / denom


def fmt(value, sig=6):
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        if not math.isfinite(float(value)):
            return ""
        return f"{float(value):.{sig}g}"
    except Exception:
        return str(value)


def decode_bytes(data):
    for encoding in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1", errors="replace")


def parse_filename(name):
    base = Path(name).name
    stem = base[:-4] if base.lower().endswith(".dat") else Path(base).stem
    m = re.search(r"IV_(positive|negative)_(.+?)(?:__\d+)?$", stem, flags=re.IGNORECASE)
    if m:
        direction = m.group(1).lower()
        raw = m.group(2)
    else:
        direction = "unknown"
        raw = stem
    condition = ""
    contact_id = raw
    # Keep condition separate so positive/negative pairs are grouped correctly.
    for cond in ("dark", "light", "illum", "illuminated"):
        token = "_" + cond
        if token in contact_id:
            condition = cond
            contact_id = contact_id.replace(token, "")
    contact_id = contact_id.strip("_") or raw
    return direction, contact_id, condition


def parse_dat_text(text, source_name):
    metadata = {}
    points = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        # Numeric data are normally: voltage/current separated by tab or spaces.
        parts = re.split(r"[\t; ,]+", line)
        if len(parts) >= 2:
            v = safe_float(parts[0])
            i = safe_float(parts[1])
            if v is not None and i is not None:
                points.append((v, i))
                continue

        # Header metadata are normally: key<TAB>value.
        if "\t" in raw_line:
            key, val = raw_line.split("\t", 1)
            metadata[key.strip()] = val.strip()
        else:
            metadata[line] = ""
    return {"name": Path(source_name).name, "metadata": metadata, "points": points}


def read_input_files(input_path):
    """Return list of parsed .dat file records from a folder or zip file."""
    input_path = Path(input_path)
    records = []

    if input_path.is_file() and input_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(input_path, "r") as zf:
            for name in sorted(zf.namelist()):
                base = Path(name).name
                if not base.lower().endswith(".dat"):
                    continue
                if base.startswith("_") or base.lower() == "_report.txt":
                    continue
                text = decode_bytes(zf.read(name))
                records.append(parse_dat_text(text, base))
        return records

    if input_path.is_dir():
        pattern = "**/*.dat" if RECURSIVE_SEARCH else "*.dat"
        for file_path in sorted(input_path.glob(pattern)):
            base = file_path.name
            if base.startswith("_") or base.lower() == "_report.txt":
                continue
            text = decode_bytes(file_path.read_bytes())
            records.append(parse_dat_text(text, str(file_path)))
        return records

    raise FileNotFoundError(f"Input path not found or unsupported: {input_path}")


def select_path_with_gui_or_prompt():
    """Ask the user to select a folder or zip if the configured path is not found."""
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox

        root = tk.Tk()
        root.withdraw()
        messagebox.showinfo(
            "Select IV data",
            "The configured DATA_PATH was not found. Select your IV data folder or IV.zip file.",
        )
        selected_file = filedialog.askopenfilename(
            title="Select IV.zip or cancel to choose a folder",
            filetypes=[("ZIP files", "*.zip"), ("All files", "*.*")],
        )
        if selected_file:
            return Path(selected_file)
        selected_dir = filedialog.askdirectory(title="Select IV data folder")
        if selected_dir:
            return Path(selected_dir)
    except Exception:
        pass

    print("DATA_PATH was not found.")
    entered = input("Paste the full folder path or ZIP path and press ENTER: ").strip().strip('"')
    return Path(entered)


def resolve_input_path():
    # Environment variable is useful for testing or temporary overrides.
    env_path = os.environ.get("IVR_DATA_PATH")
    candidates = []
    if env_path:
        candidates.append(Path(env_path))
    candidates.append(Path(DATA_PATH))

    script_dir = Path(__file__).resolve().parent
    candidates.extend([
        script_dir / "IV.zip",
        script_dir / "G5",
        script_dir,
    ])

    for candidate in candidates:
        try:
            if candidate.is_file() and candidate.suffix.lower() == ".zip":
                return candidate
            if candidate.is_dir():
                pattern = "**/*.dat" if RECURSIVE_SEARCH else "*.dat"
                if any(candidate.glob(pattern)):
                    return candidate
        except Exception:
            continue

    return select_path_with_gui_or_prompt()


def split_sweep_branches(points):
    """Split a voltage sweep into branches when the voltage step changes sign."""
    if len(points) < 3:
        return [points]
    branches = []
    start = 0
    previous_sign = 0
    for idx in range(1, len(points)):
        dv = points[idx][0] - points[idx - 1][0]
        sign = 1 if dv > 0 else (-1 if dv < 0 else 0)
        if sign != 0 and previous_sign != 0 and sign != previous_sign:
            branches.append(points[start:idx])
            start = idx - 1  # include turning point in the next branch
        if sign != 0:
            previous_sign = sign
    branches.append(points[start:])
    return [b for b in branches if len(b) >= 2]


def linear_regression(points):
    """Fit I = slope*V + intercept. Return slope, intercept, R^2."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    n = len(points)
    if n < 2:
        return None
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    sxx = sum((x - x_mean) ** 2 for x in xs)
    if sxx == 0:
        return None
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    fitted = [slope * x + intercept for x in xs]
    ss_res = sum((y - yhat) ** 2 for y, yhat in zip(ys, fitted))
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else 1.0
    return slope, intercept, r2


def manual_fit_range_for(file_name, contact_id):
    for key in (file_name, contact_id, "DEFAULT"):
        value = FIT_ABS_V_RANGE_BY_CONTACT.get(key)
        if value is not None:
            return value
    return None


def choose_linear_region(points, file_name, contact_id):
    """Choose points for the ohmic/linear fit and return fit information."""
    if len(points) < 2:
        raise ValueError("Not enough points for linear fit")

    manual_range = manual_fit_range_for(file_name, contact_id)
    if manual_range is not None:
        vmin_abs, vmax_abs = manual_range
        selected = [p for p in points if vmin_abs <= abs(p[0]) <= vmax_abs]
        if len(selected) < 2:
            raise ValueError(
                f"Manual fit range {manual_range} selected fewer than 2 points for {file_name}"
            )
        fit = linear_regression(selected)
        if fit is None:
            raise ValueError("Manual linear fit failed")
        slope, intercept, r2 = fit
        return {
            "selected_points": selected,
            "slope": slope,
            "intercept": intercept,
            "r2": r2,
            "fit_method": f"manual_absV_{vmin_abs:g}_to_{vmax_abs:g}",
        }

    # Automatic search over consecutive voltage windows.
    pts = sorted(points, key=lambda p: p[0])
    n_total = len(pts)
    min_points = min(max(3, AUTO_MIN_POINTS), n_total)
    max_span = abs(pts[-1][0] - pts[0][0])
    min_span = max_span * AUTO_MIN_SPAN_FRACTION

    best_good = None
    best_any = None

    for start in range(n_total):
        for end in range(start + min_points, n_total + 1):
            window = pts[start:end]
            span = abs(window[-1][0] - window[0][0])
            if max_span > 0 and span < min_span:
                continue
            fit = linear_regression(window)
            if fit is None:
                continue
            slope, intercept, r2 = fit
            if slope == 0 or not math.isfinite(slope):
                continue
            # Best available window prioritizes R^2 first.
            any_score = (r2, span / max_span if max_span else 1.0, len(window) / n_total)
            if best_any is None or any_score > best_any[0]:
                best_any = (any_score, window, slope, intercept, r2)
            # Good windows prioritize largest voltage span, then R^2.
            if r2 >= AUTO_MIN_R2:
                good_score = (span / max_span if max_span else 1.0, r2, len(window) / n_total)
                if best_good is None or good_score > best_good[0]:
                    best_good = (good_score, window, slope, intercept, r2)

    chosen = best_good if best_good is not None else best_any
    if chosen is None:
        # Last fallback: fit all branch points.
        fit = linear_regression(pts)
        if fit is None:
            raise ValueError("Automatic linear fit failed")
        slope, intercept, r2 = fit
        selected = pts
        method = "auto_all_points_fallback"
    else:
        _, selected, slope, intercept, r2 = chosen
        method = "auto_longest_R2_window" if best_good is not None else "auto_best_available_low_R2"

    return {
        "selected_points": selected,
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "fit_method": method,
    }


def get_geometry(contact_id, file_name):
    for key in (file_name, contact_id, "DEFAULT"):
        values = GEOMETRY_BY_CONTACT.get(key)
        if values:
            area = values.get("area_cm2")
            distance = values.get("distance_cm")
            if area is not None or distance is not None:
                return area, distance, key
    return None, None, ""


def literature_check(rho_ohm_cm):
    if rho_ohm_cm is None:
        return "geometry missing - rho not calculated"
    ranges = [r for r in LITERATURE_RANGES_OHM_CM if r.get("min") is not None and r.get("max") is not None]
    if not ranges:
        return "no literature ranges configured"
    global_min = min(min(r["min"], r["max"]) for r in ranges)
    global_max = max(max(r["min"], r["max"]) for r in ranges)
    matching = []
    for r in ranges:
        lo = min(r["min"], r["max"])
        hi = max(r["min"], r["max"])
        if lo <= rho_ohm_cm <= hi:
            matching.append(r["label"])
    if matching:
        return "within configured literature range: " + "; ".join(matching)
    if rho_ohm_cm < global_min:
        factor = global_min / rho_ohm_cm if rho_ohm_cm > 0 else float("inf")
        return f"below configured literature window by factor {factor:.3g}"
    factor = rho_ohm_cm / global_max
    return f"above configured literature window by factor {factor:.3g}"


def current_ratio_common_abs_voltage(pos_points, neg_points):
    """Median |I(+V)| / |I(-V)| for common absolute voltages."""
    def map_abs_v_to_abs_i(points):
        buckets = {}
        for v, current in points:
            key = round(abs(v), 9)
            if key == 0:
                continue
            buckets.setdefault(key, []).append(abs(current))
        return {k: median(vals) for k, vals in buckets.items()}

    pos = map_abs_v_to_abs_i(pos_points)
    neg = map_abs_v_to_abs_i(neg_points)
    common = sorted(set(pos).intersection(set(neg)))
    ratios = []
    for key in common:
        if neg[key] and neg[key] > 0 and pos[key] is not None:
            ratios.append(pos[key] / neg[key])
    if not ratios:
        return None, 0
    return median(ratios), len(ratios)


def discrepancy_explanation(ratio_r, current_ratio, low_r2, rho_statuses):
    comments = []
    if ratio_r is not None and ratio_r >= DIRECTION_RATIO_WARN:
        comments.append(
            "large positive/negative resistance asymmetry; possible rectifying/asymmetric contacts, "
            "contact polarity effect, ion migration/hysteresis, or non-ohmic transport"
        )
    elif ratio_r is not None and ratio_r >= 1.25:
        comments.append("moderate positive/negative resistance asymmetry")
    else:
        comments.append("positive/negative resistance values are relatively consistent")

    if current_ratio is not None:
        if current_ratio >= DIRECTION_RATIO_WARN or current_ratio <= 1.0 / DIRECTION_RATIO_WARN:
            comments.append("current at common |V| also shows voltage-direction asymmetry")

    if low_r2:
        comments.append("one or more linear fits have low R^2; inspect/adjust FIT_ABS_V_RANGE_BY_CONTACT")

    joined_status = " ".join([s for s in rho_statuses if s])
    if "above configured" in joined_status or "below configured" in joined_status:
        comments.append(
            "resistivity is outside configured literature window; check geometry, contact area, dark conditions, "
            "surface leakage, contact resistance, and whether the fitted region is truly ohmic"
        )
    return "; ".join(comments)


def analyze_record(record):
    file_name = record["name"]
    points = record["points"]
    direction, contact_id, condition = parse_filename(file_name)
    if not points:
        raise ValueError(f"No numeric IV points found in {file_name}")

    area_cm2, distance_cm, geometry_key = get_geometry(contact_id, file_name)
    branches = split_sweep_branches(points)
    results = []

    for branch_number, branch_points in enumerate(branches, start=1):
        branch_start_v = branch_points[0][0]
        branch_end_v = branch_points[-1][0]
        branch_direction = "increasing" if branch_end_v > branch_start_v else "decreasing"
        fit = choose_linear_region(branch_points, file_name, contact_id)
        selected = fit["selected_points"]
        slope = fit["slope"]
        intercept = fit["intercept"]
        r2 = fit["r2"]
        resistance_signed = 1.0 / slope if slope != 0 else None
        resistance_ohm = abs(resistance_signed) if resistance_signed is not None else None

        rho_ohm_cm = None
        if resistance_ohm is not None and area_cm2 is not None and distance_cm is not None:
            if area_cm2 > 0 and distance_cm > 0:
                rho_ohm_cm = resistance_ohm * area_cm2 / distance_cm

        warnings = []
        if r2 < POOR_FIT_R2_WARN:
            warnings.append("low_R2_fit")
        if fit["fit_method"] == "auto_best_available_low_R2":
            warnings.append("auto_fit_below_R2_threshold")
        if slope < 0:
            warnings.append("negative_dIdV_slope")
        if area_cm2 is None or distance_cm is None:
            warnings.append("geometry_missing")

        results.append({
            "file": file_name,
            "direction": direction,
            "contact_id": contact_id,
            "condition": condition,
            "branch": branch_number,
            "branch_direction": branch_direction,
            "branch_v_start_V": branch_start_v,
            "branch_v_end_V": branch_end_v,
            "fit_method": fit["fit_method"],
            "fit_v_min_V": min(p[0] for p in selected),
            "fit_v_max_V": max(p[0] for p in selected),
            "fit_abs_v_min_V": min(abs(p[0]) for p in selected),
            "fit_abs_v_max_V": max(abs(p[0]) for p in selected),
            "n_points_total_branch": len(branch_points),
            "n_points_fit": len(selected),
            "slope_A_per_V": slope,
            "intercept_A": intercept,
            "R2": r2,
            "resistance_signed_ohm": resistance_signed,
            "resistance_ohm": resistance_ohm,
            "area_cm2": area_cm2,
            "distance_cm": distance_cm,
            "geometry_key": geometry_key,
            "rho_ohm_cm": rho_ohm_cm,
            "literature_check": literature_check(rho_ohm_cm),
            "warning": ";".join(warnings),
        })
    return results


def summarize_file(branch_rows):
    r_values = [row["resistance_ohm"] for row in branch_rows if row.get("resistance_ohm") is not None]
    rho_values = [row["rho_ohm_cm"] for row in branch_rows if row.get("rho_ohm_cm") is not None]
    r_median = median(r_values)
    rho_median = median(rho_values)
    hysteresis_ratio = None
    if len(r_values) >= 2 and min(r_values) > 0:
        hysteresis_ratio = max(r_values) / min(r_values)

    first = branch_rows[0]
    warnings = []
    for row in branch_rows:
        if row.get("warning"):
            warnings.extend(row["warning"].split(";"))
    if hysteresis_ratio is not None and hysteresis_ratio >= HYSTERESIS_RATIO_WARN:
        warnings.append("branch_hysteresis_large")

    return {
        "file": first["file"],
        "direction": first["direction"],
        "contact_id": first["contact_id"],
        "condition": first["condition"],
        "n_branches": len(branch_rows),
        "median_resistance_ohm": r_median,
        "median_rho_ohm_cm": rho_median,
        "hysteresis_R_ratio_max_min": hysteresis_ratio,
        "best_branch_R2_min": min(row["R2"] for row in branch_rows),
        "literature_check": literature_check(rho_median),
        "warning": ";".join(sorted(set([w for w in warnings if w]))),
    }


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k), sig=10) for k in fieldnames})


def make_plots(records_by_name, branch_results_by_file, out_dir):
    if not MAKE_PLOTS:
        return "plotting disabled"
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        return f"matplotlib not available, skipped plots: {exc}"

    plot_dir = out_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    made = 0

    for file_name, record in records_by_name.items():
        points = record["points"]
        rows = branch_results_by_file[file_name]
        if not points:
            continue
        plt.figure(figsize=(7, 5))
        plt.plot([p[0] for p in points], [p[1] for p in points], "o-", markersize=3, label="data")
        for row in rows:
            x1 = row["fit_v_min_V"]
            x2 = row["fit_v_max_V"]
            y1 = row["slope_A_per_V"] * x1 + row["intercept_A"]
            y2 = row["slope_A_per_V"] * x2 + row["intercept_A"]
            label = f"fit branch {row['branch']}, R2={row['R2']:.3f}"
            plt.plot([x1, x2], [y1, y2], "--", linewidth=2, label=label)
        plt.xlabel("Voltage V")
        plt.ylabel("Current A")
        plt.title(file_name)
        plt.grid(True)
        plt.legend(fontsize=8)
        plt.tight_layout()
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(file_name).stem) + ".png"
        plt.savefig(plot_dir / safe_name, dpi=200)
        plt.close()
        made += 1
    return f"created {made} plots in {plot_dir}"


def build_pair_comparisons(file_summaries, records_by_name):
    # Group by contact and condition.
    groups = {}
    for summary in file_summaries:
        key = (summary["contact_id"], summary["condition"])
        groups.setdefault(key, {})[summary["direction"]] = summary

    pair_rows = []
    for (contact_id, condition), by_direction in sorted(groups.items()):
        pos = by_direction.get("positive")
        neg = by_direction.get("negative")
        if not pos or not neg:
            continue

        pos_points = records_by_name[pos["file"]]["points"]
        neg_points = records_by_name[neg["file"]]["points"]
        current_ratio, n_common = current_ratio_common_abs_voltage(pos_points, neg_points)
        r_ratio = ratio_max_min(pos["median_resistance_ohm"], neg["median_resistance_ohm"])
        rho_ratio = ratio_max_min(pos["median_rho_ohm_cm"], neg["median_rho_ohm_cm"])
        low_r2 = (
            (pos["best_branch_R2_min"] is not None and pos["best_branch_R2_min"] < POOR_FIT_R2_WARN)
            or (neg["best_branch_R2_min"] is not None and neg["best_branch_R2_min"] < POOR_FIT_R2_WARN)
        )
        explanation = discrepancy_explanation(
            r_ratio,
            current_ratio,
            low_r2,
            [pos.get("literature_check", ""), neg.get("literature_check", "")],
        )

        pair_rows.append({
            "contact_id": contact_id,
            "condition": condition,
            "positive_file": pos["file"],
            "negative_file": neg["file"],
            "R_positive_ohm": pos["median_resistance_ohm"],
            "R_negative_ohm": neg["median_resistance_ohm"],
            "R_ratio_max_min": r_ratio,
            "R_percent_difference_pos_vs_neg": percent_difference(pos["median_resistance_ohm"], neg["median_resistance_ohm"]),
            "rho_positive_ohm_cm": pos["median_rho_ohm_cm"],
            "rho_negative_ohm_cm": neg["median_rho_ohm_cm"],
            "rho_ratio_max_min": rho_ratio,
            "rho_percent_difference_pos_vs_neg": percent_difference(pos["median_rho_ohm_cm"], neg["median_rho_ohm_cm"]),
            "median_abs_current_ratio_pos_over_neg_common_absV": current_ratio,
            "n_common_abs_voltage_points": n_common,
            "explanation": explanation,
        })
    return pair_rows


def write_report(out_path, input_path, branch_rows, file_summary_rows, pair_rows, plot_status):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("IVR analysis report\n")
        f.write("===================\n\n")
        f.write(f"Input path: {input_path}\n")
        f.write(f"Files analyzed: {len(file_summary_rows)}\n")
        f.write(f"Total branch fits: {len(branch_rows)}\n")
        f.write(f"Plot status: {plot_status}\n\n")

        f.write("Method\n")
        f.write("------\n")
        f.write("Each .dat file is parsed as Voltage (V) and Current (A). If a file contains an up/down sweep, ")
        f.write("it is split into separate voltage branches. For each branch the script fits I = m*V + b ")
        f.write("in the selected linear region. Resistance is R = 1/m, reported as |R| in Ohm. ")
        f.write("Resistivity is rho = R*A/L in Ohm*cm when contact area A (cm^2) and distance/thickness L (cm) are entered.\n\n")

        f.write("Configured literature ranges, Ohm*cm\n")
        f.write("-------------------------------------\n")
        for r in LITERATURE_RANGES_OHM_CM:
            f.write(f"- {r['label']}: {r['min']:.3g} to {r['max']:.3g} Ohm*cm\n")
        f.write("\n")

        if any(row.get("rho_ohm_cm") is None for row in branch_rows):
            f.write("Geometry note\n")
            f.write("-------------\n")
            f.write("At least one result has no rho because area_cm2 and/or distance_cm is still None. ")
            f.write("Enter the measured geometry in GEOMETRY_BY_CONTACT and run the script again.\n\n")

        f.write("Per-file summary\n")
        f.write("----------------\n")
        for row in file_summary_rows:
            f.write(
                f"{row['file']}: R_median={fmt(row['median_resistance_ohm'])} Ohm, "
                f"rho_median={fmt(row['median_rho_ohm_cm'])} Ohm*cm, "
                f"hysteresis_ratio={fmt(row['hysteresis_R_ratio_max_min'])}, "
                f"R2_min={fmt(row['best_branch_R2_min'])}, "
                f"check={row['literature_check']}, warning={row['warning']}\n"
            )
        f.write("\n")

        f.write("Positive vs negative comparison\n")
        f.write("-------------------------------\n")
        if not pair_rows:
            f.write("No positive/negative pairs found.\n")
        for row in pair_rows:
            f.write(
                f"{row['contact_id']} ({row['condition'] or 'no condition'}): "
                f"R+={fmt(row['R_positive_ohm'])} Ohm, R-={fmt(row['R_negative_ohm'])} Ohm, "
                f"R ratio={fmt(row['R_ratio_max_min'])}, "
                f"R %diff={fmt(row['R_percent_difference_pos_vs_neg'])} %, "
                f"current ratio |I+|/|I-|={fmt(row['median_abs_current_ratio_pos_over_neg_common_absV'])}. "
                f"Explanation: {row['explanation']}\n"
            )
        f.write("\n")

        f.write("General discrepancy checklist\n")
        f.write("-----------------------------\n")
        f.write("If results differ from literature or differ strongly between +V and -V, check: actual area/thickness, ")
        f.write("whether the fitted region is ohmic, contact resistance, asymmetric/rectifying contacts, surface leakage, ")
        f.write("ionic migration/hysteresis in CsPbBr3, measurement noise at low current, dark/illumination conditions, ")
        f.write("and sample quality/defect density.\n")


def main():
    input_path = resolve_input_path()
    records = read_input_files(input_path)
    if not records:
        raise RuntimeError(f"No .dat files found in {input_path}")

    if Path(input_path).is_dir():
        out_dir = Path(input_path) / OUTPUT_FOLDER_NAME
    else:
        out_dir = Path(input_path).resolve().parent / OUTPUT_FOLDER_NAME
    out_dir.mkdir(parents=True, exist_ok=True)

    records_by_name = {record["name"]: record for record in records}
    branch_results_by_file = {}
    all_branch_rows = []
    file_summary_rows = []

    for record in records:
        rows = analyze_record(record)
        branch_results_by_file[record["name"]] = rows
        all_branch_rows.extend(rows)
        file_summary_rows.append(summarize_file(rows))

    pair_rows = build_pair_comparisons(file_summary_rows, records_by_name)

    branch_fields = [
        "file", "direction", "contact_id", "condition", "branch", "branch_direction",
        "branch_v_start_V", "branch_v_end_V", "fit_method", "fit_v_min_V", "fit_v_max_V",
        "fit_abs_v_min_V", "fit_abs_v_max_V", "n_points_total_branch", "n_points_fit",
        "slope_A_per_V", "intercept_A", "R2", "resistance_signed_ohm", "resistance_ohm",
        "area_cm2", "distance_cm", "geometry_key", "rho_ohm_cm", "literature_check", "warning",
    ]
    summary_fields = [
        "file", "direction", "contact_id", "condition", "n_branches", "median_resistance_ohm",
        "median_rho_ohm_cm", "hysteresis_R_ratio_max_min", "best_branch_R2_min",
        "literature_check", "warning",
    ]
    pair_fields = [
        "contact_id", "condition", "positive_file", "negative_file", "R_positive_ohm",
        "R_negative_ohm", "R_ratio_max_min", "R_percent_difference_pos_vs_neg",
        "rho_positive_ohm_cm", "rho_negative_ohm_cm", "rho_ratio_max_min",
        "rho_percent_difference_pos_vs_neg", "median_abs_current_ratio_pos_over_neg_common_absV",
        "n_common_abs_voltage_points", "explanation",
    ]

    write_csv(out_dir / "branch_fit_results.csv", all_branch_rows, branch_fields)
    write_csv(out_dir / "per_file_summary.csv", file_summary_rows, summary_fields)
    write_csv(out_dir / "positive_negative_comparison.csv", pair_rows, pair_fields)
    plot_status = make_plots(records_by_name, branch_results_by_file, out_dir)
    write_report(out_dir / "IVR_analysis_report.txt", input_path, all_branch_rows, file_summary_rows, pair_rows, plot_status)

    print("\nIVR analysis complete.")
    print(f"Input:  {input_path}")
    print(f"Output: {out_dir}")
    print("\nCreated files:")
    print("  - branch_fit_results.csv")
    print("  - per_file_summary.csv")
    print("  - positive_negative_comparison.csv")
    print("  - IVR_analysis_report.txt")
    print(f"  - {plot_status}")

    missing_geometry = any(row.get("rho_ohm_cm") is None for row in all_branch_rows)
    if missing_geometry:
        print("\nNote: resistivity rho was not calculated for every row because geometry is missing.")
        print("Edit GEOMETRY_BY_CONTACT near the top of IVR.py and run again.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\nERROR while running IVR.py:\n")
        traceback.print_exc()
    finally:
        if PAUSE_AT_END:
            try:
                input("\nPress ENTER to close...")
            except Exception:
                pass
