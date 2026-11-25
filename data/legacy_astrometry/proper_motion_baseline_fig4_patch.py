"""
DROP-IN FINAL CELL FOR proper_motion_fig3.ipynb

Run this after the existing cells have defined make_quantitative_pm_figure().
It reads the original JWST i2d primary headers, replaces the stale hard-coded
calendar-date baselines with MJD-AVG baselines, applies the same 160 pc annual-
parallax correction described in the paper, rebuilds the fake-source and
combined summaries, and finally calls the existing Figure 4 function without
changing any plotting style, colors, markers, layout, labels, or output names.
"""

from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
import math

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from IPython.display import display, Image as IPyImage


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

JULIAN_YEAR_DAYS = 365.25
PARALLAX_DISTANCE_PC = 160.0
N_COMBINED_SAMPLES_CORRECTED = 100_000
COMBINED_RNG_SEED_CORRECTED = 639188

# Define the target locally so this correction block does not depend on the
# original setup cell having been executed in the current kernel.
TARGET_RA_DEG = 150.0634341
TARGET_DEC_DEG = 2.3552216
TARGET = SkyCoord(
    ra=TARGET_RA_DEG * u.deg,
    dec=TARGET_DEC_DEG * u.deg,
    frame="icrs",
)

# The numerical correction is self-contained, but the final plotting call
# deliberately reuses the notebook's original helpers to preserve Figure 4's
# exact style. Fail early with a useful instruction after a kernel restart.
_missing_original_helpers = [
    name
    for name in ("filter_from_epoch_name", "make_quantitative_pm_figure")
    if name not in globals()
]
if _missing_original_helpers:
    raise RuntimeError(
        "Run the original proper_motion_fig3.ipynb cells first (or choose "
        "Kernel > Restart & Run All), then run this final correction cell. "
        "Missing original helper(s): " + ", ".join(_missing_original_helpers)
    )

# Optional: point this to a directory containing the original full i2d files.
# Leave as None to read only their primary headers from public MAST products.
LOCAL_ORIGINAL_I2D_DIR = globals().get("LOCAL_ORIGINAL_I2D_DIR", None)

BASELINE_PATCH_OUTDIR = Path("proper_motion_corrected_baseline")
BASELINE_PATCH_OUTDIR.mkdir(parents=True, exist_ok=True)

# These are the exact association-level products used to make the cutouts.
# MJD-AVG and JWST_X/Y/Z are retained in these original primary headers, but the
# old cutout-writing cell discarded most non-WCS metadata.
EPOCH_PRODUCTS = {
    "PRIMER_2023_F115W": "jw01837-o003_t003_nircam_clear-f115w_i2d.fits",
    "PRIMER_2023_F200W": "jw01837-o003_t003_nircam_clear-f200w_i2d.fits",
    "PRIMER_2023_F356W": "jw01837-o003_t003_nircam_clear-f356w_i2d.fits",
    "COSMOSWeb_2023_F115W": "jw01727-o051_t051_nircam_clear-f115w_i2d.fits",
    "COSMOS3D_2025_F115W": "jw05893-o014_t014_nircam_clear-f115w_i2d.fits",
    "COSMOS3D_2025_F200W": "jw05893-o014_t014_nircam_clear-f200w_i2d.fits",
    "COSMOS3D_2025_F356W": "jw05893-o014_t014_nircam_clear-f356w_i2d.fits",
}

ORIGINAL_IMAGE_SUMMARY_CANDIDATES = [
    Path("proper_motion_image_level/image_level_pm_summary.csv"),
    Path("image_level_pm_summary.csv"),
    Path("/mnt/data/image_level_pm_summary.csv"),
]

ORIGINAL_FAKE_TRIAL_CANDIDATES = [
    Path("ID639188_fake_source_PM_test/ID639188_fake_source_trials.csv"),
    Path("ID639188_fake_source_trials.csv"),
    Path("/mnt/data/ID639188_fake_source_trials.csv"),
]

CORRECTED_IMAGE_SUMMARY_PATH = (
    BASELINE_PATCH_OUTDIR / "image_level_pm_summary_corrected_baseline_parallax.csv"
)
CORRECTED_FAKE_TRIALS_PATH = (
    BASELINE_PATCH_OUTDIR / "ID639188_fake_source_trials_corrected_baseline_parallax.csv"
)
CORRECTED_FAKE_SUMMARY_PATH = (
    BASELINE_PATCH_OUTDIR / "ID639188_fake_source_recovery_summary_corrected_baseline_parallax.csv"
)
CORRECTED_COMBINED_SAMPLES_PATH = (
    BASELINE_PATCH_OUTDIR / "ID639188_fake_source_combined_pm_samples_corrected.csv"
)
CORRECTED_COMBINED_SUMMARY_PATH = (
    BASELINE_PATCH_OUTDIR / "ID639188_fake_source_combined_pm_summary_corrected.csv"
)
HEADER_AUDIT_PATH = BASELINE_PATCH_OUTDIR / "ID639188_original_i2d_header_epoch_audit.csv"
PAIR_AUDIT_PATH = BASELINE_PATCH_OUTDIR / "ID639188_corrected_pair_values.csv"


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def _first_existing_baseline_patch(paths):
    for path in paths:
        path = Path(path)
        if path.exists():
            return path
    return None


def _mast_product_url(product_filename):
    uri = f"mast:JWST/product/{product_filename}"
    return "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + quote(
        uri, safe=":/"
    )


def _read_original_primary_header(product_filename):
    """Read only the original i2d primary header, locally or by HTTP range."""
    local_path = None
    if LOCAL_ORIGINAL_I2D_DIR is not None:
        candidate = Path(LOCAL_ORIGINAL_I2D_DIR) / product_filename
        if candidate.exists():
            local_path = candidate

    if local_path is not None:
        with fits.open(local_path, memmap=False, lazy_load_hdus=True) as hdul:
            header = hdul[0].header.copy()
            if len(hdul) > 1:
                header.update(hdul[1].header)
            return header, f"local:{local_path}"

    remote_url = _mast_product_url(product_filename)
    try:
        # A JWST association i2d can be tens of GB. Request only the first
        # 256 KiB and stop at the primary-header END card; no image pixels are
        # downloaded and no fsspec/aiohttp dependency is required.
        request = Request(
            remote_url,
            headers={
                "Range": "bytes=0-262143",
                "User-Agent": "ID639188-proper-motion-header-audit/1.0",
            },
        )
        with urlopen(request, timeout=60) as response:
            first_bytes = response.read(262144)

        def parse_one_header(blob, start_byte):
            cards = [
                blob[offset : offset + 80]
                for offset in range(start_byte, len(blob), 80)
            ]
            end_card_index = next(
                index
                for index, card in enumerate(cards)
                if card[:8].decode("ascii", errors="ignore").strip() == "END"
            )
            header_bytes_used = int(
                math.ceil(((end_card_index + 1) * 80) / 2880.0) * 2880
            )
            header_text = blob[
                start_byte : start_byte + header_bytes_used
            ].decode("ascii", errors="strict")
            return fits.Header.fromstring(header_text, sep=""), header_bytes_used

        primary_header, primary_bytes = parse_one_header(first_bytes, 0)
        # JWST association i2d timing and observatory-position keywords are
        # normally in the SCI extension header, immediately after the empty
        # primary HDU. Merge it exactly as the old image loader merged headers.
        science_header, _ = parse_one_header(first_bytes, primary_bytes)
        primary_header.update(science_header)
        return primary_header, "MAST primary+SCI headers"
    except Exception as exc:
        raise RuntimeError(
            f"Could not read the primary header for {product_filename}. "
            "Set LOCAL_ORIGINAL_I2D_DIR to a folder containing the original "
            "i2d files, or check network access to MAST."
        ) from exc


def _header_midpoint_mjd(header):
    """Use MJD-AVG/EXPMID; allow explicit midpoint fallbacks with provenance."""
    for key in ("MJD-AVG", "EXPMID"):
        value = header.get(key)
        if value is not None and np.isfinite(float(value)):
            return float(value), key

    for start_key, end_key in (
        ("MJD-BEG", "MJD-END"),
        ("EXPSTART", "EXPEND"),
    ):
        start = header.get(start_key)
        end = header.get(end_key)
        if start is not None and end is not None:
            return 0.5 * (float(start) + float(end)), f"midpoint({start_key},{end_key})"

    date_beg = header.get("DATE-BEG")
    date_end = header.get("DATE-END")
    if date_beg is not None and date_end is not None:
        t1 = Time(date_beg, format="isot", scale="utc")
        t2 = Time(date_end, format="isot", scale="utc")
        return float(0.5 * (t1.mjd + t2.mjd)), "midpoint(DATE-BEG,DATE-END)"

    raise KeyError(
        "No MJD-AVG, EXPMID, MJD-BEG/MJD-END, EXPSTART/EXPEND, "
        "or DATE-BEG/DATE-END keywords were found."
    )


def _robust_sigma_baseline_patch(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    median = np.nanmedian(values)
    return float(1.4826 * np.nanmedian(np.abs(values - median)))


def _percentile_or_nan_baseline_patch(values, q):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.nanpercentile(values, q)) if values.size else np.nan


def _proper_motion_pa_vectorized(mu_ra, mu_dec):
    return (np.degrees(np.arctan2(mu_ra, mu_dec)) + 360.0) % 360.0


# -----------------------------------------------------------------------------
# 1. Audit original headers and derive the correct epochs
# -----------------------------------------------------------------------------

header_rows = []
for epoch_name, product_filename in EPOCH_PRODUCTS.items():
    header, access_source = _read_original_primary_header(product_filename)
    midpoint_mjd, midpoint_keyword = _header_midpoint_mjd(header)

    mjd_beg = float(header.get("MJD-BEG", header.get("EXPSTART", np.nan)))
    mjd_end = float(header.get("MJD-END", header.get("EXPEND", np.nan)))
    if np.isfinite(mjd_beg) and np.isfinite(mjd_end):
        if not (mjd_beg <= midpoint_mjd <= mjd_end):
            raise ValueError(
                f"Midpoint for {epoch_name} is outside the header exposure range."
            )

    xyz = [header.get(key) for key in ("JWST_X", "JWST_Y", "JWST_Z")]
    if any(value is None for value in xyz):
        raise KeyError(
            f"JWST_X/Y/Z are missing from {product_filename}; the parallax "
            "correction described in the paper cannot be reproduced."
        )

    header_rows.append(
        {
            "epoch": epoch_name,
            "product": product_filename,
            "midpoint_keyword": midpoint_keyword,
            "mjd_avg": midpoint_mjd,
            "midpoint_utc": Time(midpoint_mjd, format="mjd", scale="utc").isot,
            "mjd_beg": mjd_beg,
            "mjd_end": mjd_end,
            "jwst_x_km": float(xyz[0]),
            "jwst_y_km": float(xyz[1]),
            "jwst_z_km": float(xyz[2]),
            "header_access": access_source,
        }
    )

header_audit = pd.DataFrame(header_rows).sort_values("mjd_avg").reset_index(drop=True)
header_audit.to_csv(HEADER_AUDIT_PATH, index=False)
epoch_metadata = header_audit.set_index("epoch")

if set(EPOCH_PRODUCTS) != set(epoch_metadata.index):
    raise RuntimeError("The header audit did not return exactly one row per epoch.")

print("Original i2d header midpoint audit:")
print(
    header_audit[
        ["epoch", "product", "midpoint_keyword", "mjd_avg", "midpoint_utc"]
    ].to_string(index=False)
)


# -----------------------------------------------------------------------------
# 2. Recompute baselines, parallax-corrected displacements, and proper motions
# -----------------------------------------------------------------------------

original_image_summary_path = _first_existing_baseline_patch(
    ORIGINAL_IMAGE_SUMMARY_CANDIDATES
)
if original_image_summary_path is None:
    raise FileNotFoundError(
        "Could not find the original image_level_pm_summary.csv. Run the image-"
        "level proper-motion notebook first."
    )

original_image_summary = pd.read_csv(original_image_summary_path)
required_image_columns = {
    "epoch1",
    "epoch2",
    "cand_corr_dra_mas",
    "cand_corr_ddec_mas",
}
missing_image_columns = required_image_columns - set(original_image_summary.columns)
if missing_image_columns:
    raise KeyError(f"Missing image-summary columns: {sorted(missing_image_columns)}")

corrected_image_summary = original_image_summary.copy()
corrected_image_summary["dt_yr_old"] = pd.to_numeric(
    corrected_image_summary.get("dt_yr", np.nan), errors="coerce"
)

target_ra_rad = float(TARGET.ra.radian)
target_dec_rad = float(TARGET.dec.radian)
au_km = float(u.au.to(u.km))
parallax_mas = 1000.0 / PARALLAX_DISTANCE_PC


def _parallax_factors_from_epoch(epoch_name):
    row = epoch_metadata.loc[epoch_name]
    x = float(row["jwst_x_km"]) / au_km
    y = float(row["jwst_y_km"]) / au_km
    z = float(row["jwst_z_km"]) / au_km
    p_ra = x * math.sin(target_ra_rad) - y * math.cos(target_ra_rad)
    p_dec = (
        x * math.cos(target_ra_rad) * math.sin(target_dec_rad)
        + y * math.sin(target_ra_rad) * math.sin(target_dec_rad)
        - z * math.cos(target_dec_rad)
    )
    return float(p_ra), float(p_dec)


for idx, row in corrected_image_summary.iterrows():
    epoch1 = str(row["epoch1"])
    epoch2 = str(row["epoch2"])
    if epoch1 not in epoch_metadata.index or epoch2 not in epoch_metadata.index:
        raise KeyError(f"No original-header metadata for {epoch1} -> {epoch2}")

    mjd1 = float(epoch_metadata.loc[epoch1, "mjd_avg"])
    mjd2 = float(epoch_metadata.loc[epoch2, "mjd_avg"])
    baseline_days = mjd2 - mjd1
    dt_yr = baseline_days / JULIAN_YEAR_DAYS
    if not np.isfinite(dt_yr) or dt_yr <= 0:
        raise ValueError(f"Non-positive header baseline for {epoch1} -> {epoch2}")

    p_ra1, p_dec1 = _parallax_factors_from_epoch(epoch1)
    p_ra2, p_dec2 = _parallax_factors_from_epoch(epoch2)
    parallax_dra = parallax_mas * (p_ra2 - p_ra1)
    parallax_ddec = parallax_mas * (p_dec2 - p_dec1)

    local_dra = float(row["cand_corr_dra_mas"])
    local_ddec = float(row["cand_corr_ddec_mas"])
    intrinsic_dra = local_dra - parallax_dra
    intrinsic_ddec = local_ddec - parallax_ddec

    corrected_image_summary.loc[idx, "mjd_avg_epoch1"] = mjd1
    corrected_image_summary.loc[idx, "mjd_avg_epoch2"] = mjd2
    corrected_image_summary.loc[idx, "baseline_days"] = baseline_days
    corrected_image_summary.loc[idx, "dt_yr"] = dt_yr
    corrected_image_summary.loc[idx, "parallax_distance_pc"] = PARALLAX_DISTANCE_PC
    corrected_image_summary.loc[idx, "parallax_dra_mas"] = parallax_dra
    corrected_image_summary.loc[idx, "parallax_ddec_mas"] = parallax_ddec
    corrected_image_summary.loc[idx, "cand_localframe_dra_mas"] = local_dra
    corrected_image_summary.loc[idx, "cand_localframe_ddec_mas"] = local_ddec
    corrected_image_summary.loc[idx, "cand_corr_dra_mas"] = intrinsic_dra
    corrected_image_summary.loc[idx, "cand_corr_ddec_mas"] = intrinsic_ddec
    corrected_image_summary.loc[idx, "cand_total_mas"] = math.hypot(
        intrinsic_dra, intrinsic_ddec
    )
    corrected_image_summary.loc[idx, "mu_ra_masyr"] = intrinsic_dra / dt_yr
    corrected_image_summary.loc[idx, "mu_dec_masyr"] = intrinsic_ddec / dt_yr
    corrected_image_summary.loc[idx, "mu_total_masyr"] = (
        math.hypot(intrinsic_dra, intrinsic_ddec) / dt_yr
    )

corrected_image_summary.to_csv(CORRECTED_IMAGE_SUMMARY_PATH, index=False)

pair_audit_columns = [
    "epoch1",
    "epoch2",
    "dt_yr_old",
    "dt_yr",
    "baseline_days",
    "parallax_dra_mas",
    "parallax_ddec_mas",
    "mu_ra_masyr",
    "mu_dec_masyr",
    "mu_total_masyr",
]
pair_audit = corrected_image_summary[pair_audit_columns].copy()
pair_audit.to_csv(PAIR_AUDIT_PATH, index=False)

print("\nCorrected pair values (MJD-AVG baseline; d = 160 pc parallax):")
print(pair_audit.to_string(index=False))


# -----------------------------------------------------------------------------
# 3. Rebuild fake-source errors in mas/yr using the corrected pair baselines
# -----------------------------------------------------------------------------

original_fake_trials_path = _first_existing_baseline_patch(
    ORIGINAL_FAKE_TRIAL_CANDIDATES
)
if original_fake_trials_path is None:
    raise FileNotFoundError(
        "The raw fake-source trial CSV is required for an exact corrected Figure 4. "
        "Run fake_source_proper_motion.ipynb first; do not reuse only the old "
        "mas/yr summary because its errors also contain the stale baseline."
    )

corrected_trials = pd.read_csv(original_fake_trials_path)
required_trial_columns = {
    "epoch1",
    "epoch2",
    "mode",
    "input_corr_dra_mas",
    "input_corr_ddec_mas",
    "rec_corr_dra_mas",
    "rec_corr_ddec_mas",
}
missing_trial_columns = required_trial_columns - set(corrected_trials.columns)
if missing_trial_columns:
    raise KeyError(f"Missing fake-trial columns: {sorted(missing_trial_columns)}")

corrected_trials["dt_yr_old"] = pd.to_numeric(
    corrected_trials.get("dt_yr", np.nan), errors="coerce"
)
corrected_trials["input_localframe_dra_mas"] = corrected_trials[
    "input_corr_dra_mas"
]
corrected_trials["input_localframe_ddec_mas"] = corrected_trials[
    "input_corr_ddec_mas"
]
corrected_trials["rec_localframe_dra_mas"] = corrected_trials[
    "rec_corr_dra_mas"
]
corrected_trials["rec_localframe_ddec_mas"] = corrected_trials[
    "rec_corr_ddec_mas"
]

for _, pair_row in corrected_image_summary.iterrows():
    epoch1 = str(pair_row["epoch1"])
    epoch2 = str(pair_row["epoch2"])
    pair_mask = (
        corrected_trials["epoch1"].astype(str).eq(epoch1)
        & corrected_trials["epoch2"].astype(str).eq(epoch2)
    )
    if not pair_mask.any():
        raise ValueError(f"No fake trials found for {epoch1} -> {epoch2}")

    dt_yr = float(pair_row["dt_yr"])
    parallax_dra = float(pair_row["parallax_dra_mas"])
    parallax_ddec = float(pair_row["parallax_ddec_mas"])
    moving_mask = pair_mask & corrected_trials["mode"].astype(str).str.lower().eq(
        "moving"
    )

    # Only moving fakes emulate the target's apparent displacement. The
    # stationary injections remain a detector/registration null test at zero.
    corrected_trials.loc[moving_mask, "input_corr_dra_mas"] -= parallax_dra
    corrected_trials.loc[moving_mask, "input_corr_ddec_mas"] -= parallax_ddec
    corrected_trials.loc[moving_mask, "rec_corr_dra_mas"] -= parallax_dra
    corrected_trials.loc[moving_mask, "rec_corr_ddec_mas"] -= parallax_ddec

    corrected_trials.loc[pair_mask, "dt_yr"] = dt_yr

    input_dra = corrected_trials.loc[pair_mask, "input_corr_dra_mas"].to_numpy(float)
    input_ddec = corrected_trials.loc[pair_mask, "input_corr_ddec_mas"].to_numpy(float)
    rec_dra = corrected_trials.loc[pair_mask, "rec_corr_dra_mas"].to_numpy(float)
    rec_ddec = corrected_trials.loc[pair_mask, "rec_corr_ddec_mas"].to_numpy(float)

    err_dra = rec_dra - input_dra
    err_ddec = rec_ddec - input_ddec
    input_mu_ra = input_dra / dt_yr
    input_mu_dec = input_ddec / dt_yr
    rec_mu_ra = rec_dra / dt_yr
    rec_mu_dec = rec_ddec / dt_yr
    input_mu_total = np.hypot(input_mu_ra, input_mu_dec)
    rec_mu_total = np.hypot(rec_mu_ra, rec_mu_dec)

    corrected_trials.loc[pair_mask, "input_mu_ra_masyr"] = input_mu_ra
    corrected_trials.loc[pair_mask, "input_mu_dec_masyr"] = input_mu_dec
    corrected_trials.loc[pair_mask, "input_mu_total_masyr"] = input_mu_total
    corrected_trials.loc[pair_mask, "rec_mu_ra_masyr"] = rec_mu_ra
    corrected_trials.loc[pair_mask, "rec_mu_dec_masyr"] = rec_mu_dec
    corrected_trials.loc[pair_mask, "rec_mu_total_masyr"] = rec_mu_total
    corrected_trials.loc[pair_mask, "rec_pa_deg_east_of_north"] = (
        _proper_motion_pa_vectorized(rec_mu_ra, rec_mu_dec)
    )
    corrected_trials.loc[pair_mask, "err_corr_dra_mas"] = err_dra
    corrected_trials.loc[pair_mask, "err_corr_ddec_mas"] = err_ddec
    corrected_trials.loc[pair_mask, "err_corr_total_mas"] = np.hypot(
        err_dra, err_ddec
    )
    corrected_trials.loc[pair_mask, "err_mu_ra_masyr"] = err_dra / dt_yr
    corrected_trials.loc[pair_mask, "err_mu_dec_masyr"] = err_ddec / dt_yr
    corrected_trials.loc[pair_mask, "err_mu_vector_masyr"] = (
        np.hypot(err_dra, err_ddec) / dt_yr
    )
    corrected_trials.loc[pair_mask, "err_mu_total_masyr"] = (
        rec_mu_total - input_mu_total
    )

corrected_trials.to_csv(CORRECTED_FAKE_TRIALS_PATH, index=False)


def _build_corrected_recovery_summary(trials, image_summary):
    summary_rows = []
    statistic_columns = [
        "err_corr_dra_mas",
        "err_corr_ddec_mas",
        "err_corr_total_mas",
        "err_mu_ra_masyr",
        "err_mu_dec_masyr",
        "err_mu_vector_masyr",
        "err_mu_total_masyr",
    ]

    for (epoch1, epoch2, mode), group in trials.groupby(
        ["epoch1", "epoch2", "mode"], sort=True
    ):
        image_match = image_summary[
            image_summary["epoch1"].astype(str).eq(str(epoch1))
            & image_summary["epoch2"].astype(str).eq(str(epoch2))
        ]
        if len(image_match) != 1:
            raise ValueError(
                f"Expected one corrected image row for {epoch1} -> {epoch2}; "
                f"found {len(image_match)}."
            )
        image_row = image_match.iloc[0]
        result = {
            "epoch1": epoch1,
            "epoch2": epoch2,
            "filter": filter_from_epoch_name(str(epoch1)),
            "mode": mode,
            "n_fake": len(group),
            "dt_yr": float(image_row["dt_yr"]),
            "observed_corr_dra_mas": float(image_row["cand_corr_dra_mas"]),
            "observed_corr_ddec_mas": float(image_row["cand_corr_ddec_mas"]),
            "observed_mu_ra_masyr": float(image_row["mu_ra_masyr"]),
            "observed_mu_dec_masyr": float(image_row["mu_dec_masyr"]),
            "observed_mu_total_masyr": float(image_row["mu_total_masyr"]),
        }

        for column in statistic_columns:
            values = pd.to_numeric(group[column], errors="coerce").to_numpy(float)
            result[f"{column}_median"] = _percentile_or_nan_baseline_patch(
                values, 50
            )
            result[f"{column}_sigma_mad"] = _robust_sigma_baseline_patch(values)
            result[f"{column}_p16"] = _percentile_or_nan_baseline_patch(values, 16)
            result[f"{column}_p84"] = _percentile_or_nan_baseline_patch(values, 84)
            result[f"{column}_p2p5"] = _percentile_or_nan_baseline_patch(
                values, 2.5
            )
            result[f"{column}_p97p5"] = _percentile_or_nan_baseline_patch(
                values, 97.5
            )

        if str(mode).lower() == "stationary":
            false_mu = pd.to_numeric(
                group["rec_mu_total_masyr"], errors="coerce"
            ).to_numpy(float)
            p95 = _percentile_or_nan_baseline_patch(false_mu, 95)
            result["stationary_false_mu_p95_masyr"] = p95
            result["observed_over_stationary_p95"] = (
                float(image_row["mu_total_masyr"]) / p95
                if np.isfinite(p95) and p95 > 0
                else np.nan
            )

        summary_rows.append(result)

    return pd.DataFrame(summary_rows)


corrected_fake_summary = _build_corrected_recovery_summary(
    corrected_trials, corrected_image_summary
)
corrected_fake_summary.to_csv(CORRECTED_FAKE_SUMMARY_PATH, index=False)


# -----------------------------------------------------------------------------
# 4. Rebuild the combined three-filter Monte Carlo result
# -----------------------------------------------------------------------------

primary_pairs_for_combination = [
    ("PRIMER_2023_F115W", "COSMOS3D_2025_F115W"),
    ("PRIMER_2023_F200W", "COSMOS3D_2025_F200W"),
    ("PRIMER_2023_F356W", "COSMOS3D_2025_F356W"),
]

combination_inputs = []
raw_weights = []
for epoch1, epoch2 in primary_pairs_for_combination:
    image_match = corrected_image_summary[
        corrected_image_summary["epoch1"].astype(str).eq(epoch1)
        & corrected_image_summary["epoch2"].astype(str).eq(epoch2)
    ]
    fake_match = corrected_trials[
        corrected_trials["epoch1"].astype(str).eq(epoch1)
        & corrected_trials["epoch2"].astype(str).eq(epoch2)
        & corrected_trials["mode"].astype(str).str.lower().eq("moving")
    ].reset_index(drop=True)

    if len(image_match) != 1 or len(fake_match) < 10:
        raise ValueError(
            f"Incomplete corrected combination inputs for {epoch1} -> {epoch2}."
        )

    image_row = image_match.iloc[0]
    # The manuscript's combined amplitude is the inverse-variance combination
    # of the three scalar |mu| measurements. Use the empirical total-amplitude
    # scatter from the moving fake recoveries for those weights. Center the
    # recovery errors below because these injections are being used to estimate
    # scatter, not to impose a centroid-bias correction on the measured motion.
    sigma_total = _robust_sigma_baseline_patch(
        fake_match["err_mu_total_masyr"]
    )
    if not np.isfinite(sigma_total) or sigma_total <= 0:
        raise ValueError(f"Invalid fake-source scatter for {epoch1} -> {epoch2}")

    raw_weights.append(1.0 / sigma_total**2)
    combination_inputs.append((image_row, fake_match))

pair_weights = np.asarray(raw_weights, dtype=float)
pair_weights /= pair_weights.sum()

rng = np.random.default_rng(COMBINED_RNG_SEED_CORRECTED)
combined_mu_ra = np.zeros(N_COMBINED_SAMPLES_CORRECTED, dtype=float)
combined_mu_dec = np.zeros(N_COMBINED_SAMPLES_CORRECTED, dtype=float)
combined_mu_total = np.zeros(N_COMBINED_SAMPLES_CORRECTED, dtype=float)

for weight, (image_row, fake_match) in zip(pair_weights, combination_inputs):
    indices = rng.integers(
        0, len(fake_match), size=N_COMBINED_SAMPLES_CORRECTED
    )
    all_error_ra = fake_match["err_mu_ra_masyr"].to_numpy(float)
    all_error_dec = fake_match["err_mu_dec_masyr"].to_numpy(float)
    all_error_total = fake_match["err_mu_total_masyr"].to_numpy(float)
    error_ra = all_error_ra[indices] - np.nanmedian(all_error_ra)
    error_dec = all_error_dec[indices] - np.nanmedian(all_error_dec)
    error_total = all_error_total[indices] - np.nanmedian(all_error_total)
    combined_mu_ra += weight * (float(image_row["mu_ra_masyr"]) - error_ra)
    combined_mu_dec += weight * (float(image_row["mu_dec_masyr"]) - error_dec)
    combined_mu_total += weight * (
        float(image_row["mu_total_masyr"]) - error_total
    )

combined_pa = _proper_motion_pa_vectorized(combined_mu_ra, combined_mu_dec)

corrected_combined_samples = pd.DataFrame(
    {
        "sample": np.arange(N_COMBINED_SAMPLES_CORRECTED),
        "mu_ra_masyr": combined_mu_ra,
        "mu_dec_masyr": combined_mu_dec,
        "mu_total_masyr": combined_mu_total,
        "pa_deg_east_of_north": combined_pa,
    }
)
corrected_combined_samples.to_csv(CORRECTED_COMBINED_SAMPLES_PATH, index=False)

combined_summary_rows = []
for column in [
    "mu_ra_masyr",
    "mu_dec_masyr",
    "mu_total_masyr",
    "pa_deg_east_of_north",
]:
    values = corrected_combined_samples[column].to_numpy(float)
    combined_summary_rows.append(
        {
            "quantity": column,
            "median": _percentile_or_nan_baseline_patch(values, 50),
            "p16": _percentile_or_nan_baseline_patch(values, 16),
            "p84": _percentile_or_nan_baseline_patch(values, 84),
            "p2p5": _percentile_or_nan_baseline_patch(values, 2.5),
            "p97p5": _percentile_or_nan_baseline_patch(values, 97.5),
            "sigma_mad": _robust_sigma_baseline_patch(values),
        }
    )

corrected_combined_summary = pd.DataFrame(combined_summary_rows)
corrected_combined_summary.to_csv(CORRECTED_COMBINED_SUMMARY_PATH, index=False)

print("\nCorrected combined main-pair result:")
print(corrected_combined_summary.to_string(index=False))


# -----------------------------------------------------------------------------
# 5. Regenerate Figure 4 using the original plotting function unchanged
# -----------------------------------------------------------------------------

# The existing helper functions read these candidate lists at call time. Point
# them at corrected tables; make_quantitative_pm_figure() itself is untouched.
SUMMARY_CSV_CANDIDATES = [CORRECTED_IMAGE_SUMMARY_PATH]
FAKE_SUMMARY_CANDIDATES = [CORRECTED_FAKE_SUMMARY_PATH]
COMBINED_SUMMARY_CANDIDATES = [CORRECTED_COMBINED_SUMMARY_PATH]
REQUIRE_FAKE_SOURCE_ERRORS = True

print("\nRegenerating Figure 4 with the original style and output filename...")
quantitative_pdf, quantitative_png = make_quantitative_pm_figure()

print("\nCorrected Figure 4 files:")
print(" ", quantitative_pdf)
print(" ", quantitative_png)
print("Header audit:", HEADER_AUDIT_PATH)
print("Pair audit:", PAIR_AUDIT_PATH)

display(IPyImage(filename=str(quantitative_png)))
