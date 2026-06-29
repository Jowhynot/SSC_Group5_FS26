from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


# ---------------- USER SETTINGS ----------------
DATA_FOLDER = Path(
    r"C:\Users\jsjan\Downloads\SSC Data&Codes\Group 5_abs_photolumin\Group 5"
)

# Use None to automatically analyze the newest CSV file in DATA_FOLDER.
# Or put a specific file path here, for example:
# CSV_FILE = DATA_FOLDER / "your_measurement.csv"
CSV_FILE = None

# Normalization options:
# "max"    -> S1c / max(S1c), common for PL spectra
# "minmax" -> (S1c - min) / (max - min), useful if there is negative background
NORMALIZATION = "max"
# ------------------------------------------------


def find_newest_csv(folder: Path) -> Path:
    """Find the newest CSV file in the measurement folder."""
    csv_files = list(folder.glob("*.csv"))

    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in: {folder}")

    return max(csv_files, key=lambda file: file.stat().st_mtime)


def load_s1c_data(csv_path: Path) -> pd.DataFrame:
    """
    Load wavelength and S1c data from a Horiba-style CSV file.

    The uploaded example has:
        Row 1: Wavelength,S1c,Wavelength,S1,Wavelength,R1
        Row 2: nm,CPS,nm,CPS,nm,MicroAmps

    This function drops the units row automatically.
    """
    raw = pd.read_csv(csv_path)

    # First wavelength column and S1c channel
    wavelength = pd.to_numeric(raw.iloc[:, 0], errors="coerce")
    s1c = pd.to_numeric(raw["S1c"], errors="coerce")

    df = pd.DataFrame({
        "Wavelength_nm": wavelength,
        "S1c_CPS": s1c
    })

    # Remove unit row and any invalid rows
    df = df.dropna().reset_index(drop=True)

    return df


def normalize_signal(df: pd.DataFrame, mode: str = "max") -> pd.DataFrame:
    """Normalize the S1c signal."""
    y = df["S1c_CPS"]

    if mode == "max":
        df["S1c_normalized"] = y / y.max()

    elif mode == "minmax":
        df["S1c_normalized"] = (y - y.min()) / (y.max() - y.min())

    else:
        raise ValueError("NORMALIZATION must be either 'max' or 'minmax'.")

    return df


def find_pl_peak(df: pd.DataFrame):
    """Find the wavelength where normalized S1c is maximum."""
    peak_index = df["S1c_normalized"].idxmax()
    peak_wavelength = df.loc[peak_index, "Wavelength_nm"]
    peak_intensity = df.loc[peak_index, "S1c_CPS"]

    return peak_wavelength, peak_intensity


def save_plot(df: pd.DataFrame, csv_path: Path, peak_wavelength: float) -> Path:
    """Create and save normalized PL plot."""
    output_plot = DATA_FOLDER / f"{csv_path.stem}_S1c_normalized_PL.png"

    plt.figure(figsize=(8, 5))
    plt.plot(
        df["Wavelength_nm"],
        df["S1c_normalized"],
        linewidth=1.8,
        label="Normalized S1c"
    )

    plt.axvline(
        peak_wavelength,
        linestyle="--",
        linewidth=1.2,
        label=f"Peak: {peak_wavelength:.1f} nm"
    )

    plt.xlabel("Wavelength / nm")
    plt.ylabel("Normalized PL intensity / -")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(output_plot, dpi=300)
    plt.show()

    return output_plot


def save_normalized_data(df: pd.DataFrame, csv_path: Path) -> Path:
    """Save normalized data as CSV."""
    output_csv = DATA_FOLDER / f"{csv_path.stem}_S1c_normalized.csv"
    df.to_csv(output_csv, index=False)
    return output_csv


def main():
    print("Starting photoluminescence data analysis...")

    if CSV_FILE is None:
        csv_path = find_newest_csv(DATA_FOLDER)
    else:
        csv_path = Path(CSV_FILE)

    print(f"Using file: {csv_path}")

    df = load_s1c_data(csv_path)
    df = normalize_signal(df, NORMALIZATION)

    peak_wavelength, peak_intensity = find_pl_peak(df)

    plot_path = save_plot(df, csv_path, peak_wavelength)
    data_path = save_normalized_data(df, csv_path)

    print()
    print("Done.")
    print(f"PL peak wavelength: {peak_wavelength:.1f} nm")
    print(f"Raw S1c intensity at peak: {peak_intensity:.3f} CPS")
    print(f"Saved plot to: {plot_path}")
    print(f"Saved normalized data to: {data_path}")


if __name__ == "__main__":
    main()