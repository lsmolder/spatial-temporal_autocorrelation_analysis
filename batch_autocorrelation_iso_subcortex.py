"""
batch_autocorrelation_iso_subcortex.py

Batch spatial and temporal autocorrelation — ISO-SUBCORTEX CROSS EDGES ONLY.

Analyses the off-diagonal block of the full 98-parcel FC matrix:
    rows = isocortical regions
    cols = subcortical regions

This captures connectivity BETWEEN iso and subcortex (e.g. thalamocortical,
cortico-limbic), which is not captured by the iso-only or subcortex-only scripts.

SA is computed using all pairwise distances between iso and sub centroids,
with FC values taken from the rectangular cross-block (no upper-triangle
extraction needed since the block is already non-redundant).

TA is computed per-region across all iso+sub regions combined, then reported
as a single global mean and per-region CSV.

Two modes:

  COHORT MODE (lecanemab project):
    python batch_autocorrelation_iso_subcortex.py \
        --parent_dir /path/to/cohorts_parent \
        --hierarchy  /path/to/ambaCCFv3_atlas_annotationsTable.csv \
        --output_dir /path/to/output

  DIRECT MODE (flat datasets, e.g. aya_project):
    python batch_autocorrelation_iso_subcortex.py \
        --timeseries_dir /path/to/timeseries_csv \
        --fc_dir         /path/to/fc_matrix \
        --atlas          /path/to/atlas.nii.gz \
        --hierarchy      /path/to/ambaCCFv3_atlas_annotationsTable.csv \
        --output_dir     /path/to/output \
        --ignore_direction
"""

import os
import argparse
import glob
import re
import numpy as np
import pandas as pd
import nibabel as nib

from calculate_autocorrelation import extract_centroids, calculate_ta
from region_filter import get_isocortex_labels, get_subcortex_labels, filter_fc_cross


# =============================================================================
# COHORT CONFIGURATION
# =============================================================================

COHORT_NAMES = [
    'cohort_100',
    'cohort_200_AP',
    'cohort_200_PA',
    'cohort_300_AP',
    'cohort_300_PA',
]

COHORTS_WITH_DIRECTION = {
    'cohort_100', 'cohort_200_AP', 'cohort_200_PA',
    'cohort_300_AP', 'cohort_300_PA'
}


# =============================================================================
# BIDS PARSING
# =============================================================================

def parse_bids(name, has_direction=True):
    sub_m = re.search(r'(sub-[a-zA-Z0-9]+)', name)
    ses_m = re.search(r'(ses-[a-zA-Z0-9]+)', name)
    run_m = re.search(r'run-([a-zA-Z0-9]+)', name)

    sub = sub_m.group(1) if sub_m else 'unknown'
    ses = ses_m.group(1) if ses_m else 'unknown'
    run = run_m.group(1) if run_m else 'unknown'

    if has_direction:
        dir_pa = re.search(r'dir-PA', name, re.IGNORECASE)
        dir_ap = re.search(r'dir-AP', name, re.IGNORECASE)
        direction = 'PA' if dir_pa else ('AP' if dir_ap else 'unknown')
    else:
        direction = 'none'

    return sub, ses, direction, run


def output_run_label(direction, run, cohort=None):
    if cohort and direction == 'PA' and cohort != 'cohort_100':
        return 'run-2'
    return f'run-{run}' if not str(run).startswith('run-') else run


# =============================================================================
# FC MATRIX FINDER
# =============================================================================

def find_fc_file(fc_root, sub, ses, direction, run):
    if not os.path.exists(fc_root):
        return None
    run_token = f'run-{run}' if not str(run).startswith('run-') else run
    for d in os.listdir(fc_root):
        dpath = os.path.join(fc_root, d)
        if not os.path.isdir(dpath):
            continue
        if sub not in d or ses not in d or run_token not in d:
            continue
        if direction not in ('unknown', 'none') and f'dir-{direction}' not in d:
            continue
        csv_files = glob.glob(os.path.join(dpath, '*fc_matrix.csv'))
        if csv_files:
            return csv_files[0]
    return None


# =============================================================================
# SORTING
# =============================================================================

def sort_summary(df):
    def ses_key(s):
        d = re.findall(r'\d+', str(s))
        return int(d[0]) if d else float('inf')
    def run_key(r):
        d = re.findall(r'\d+', str(r))
        return int(d[0]) if d else float('inf')
    df['_ses_sort'] = df['Session'].apply(ses_key)
    df['_run_sort'] = df['Run'].apply(run_key)
    df = df.sort_values(['Subject', '_ses_sort', '_run_sort'])
    return df.drop(columns=['_ses_sort', '_run_sort'])


# =============================================================================
# CROSS-BLOCK SA
# =============================================================================

def calculate_sa_cross(fc_cross, iso_centroids, sub_centroids):
    """
    Spatial autocorrelation for the rectangular iso x sub FC block.

    Distances are computed between every iso centroid and every sub centroid
    (i.e. the full cross-product of distances, shape n_iso x n_sub).
    FC values are taken directly from the cross-block — no upper triangle
    extraction needed since rows != cols so every entry is a unique pair.

    Uses 1mm bins matching the original calculate_sa methodology.

    Returns:
        sa_lambda, sa_inf, bin_centers, bin_means, fit_y
    """
    from scipy.optimize import curve_fit

    print("\nCalculating cross-block SA (iso x sub)...")

    # Pairwise distances: n_iso x n_sub
    n_iso = iso_centroids.shape[0]
    n_sub = sub_centroids.shape[0]
    dist_cross = np.zeros((n_iso, n_sub))
    for i in range(n_iso):
        for j in range(n_sub):
            dist_cross[i, j] = np.linalg.norm(
                iso_centroids[i] - sub_centroids[j])

    # Flatten both matrices — every entry is a unique iso-sub pair
    flat_dist = dist_cross.flatten()
    flat_fc   = fc_cross.flatten()

    # Remove NaNs
    valid = ~np.isnan(flat_fc) & ~np.isnan(flat_dist)
    flat_dist = flat_dist[valid]
    flat_fc   = flat_fc[valid]

    # 1mm bins
    max_dist  = int(np.ceil(flat_dist.max()))
    bins      = np.arange(0, max_dist + 1, 1)
    bin_centers, bin_means = [], []

    for i in range(len(bins) - 1):
        mask = (flat_dist >= bins[i]) & (flat_dist < bins[i + 1])
        vals = flat_fc[mask]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            bin_centers.append((bins[i] + bins[i + 1]) / 2)
            bin_means.append(float(np.mean(vals)))

    bin_centers = np.array(bin_centers)
    bin_means   = np.array(bin_means)

    print(f"  Data points for fitting: {len(bin_centers)} bins")

    if len(bin_centers) < 3:
        print("  Not enough bins to fit curve.")
        return np.nan, np.nan, bin_centers, bin_means, None

    def exponential_decay(x, sa_lambda, sa_inf):
        return sa_inf + (1 - sa_inf) * np.exp(-x / sa_lambda)

    try:
        popt, _ = curve_fit(
            exponential_decay, bin_centers, bin_means,
            p0=[10.0, 0.1],
            bounds=([0.001, -1.0], [200.0, 1.0]),
            maxfev=5000
        )
        sa_lambda, sa_inf = popt
        print(f"  SA-lambda: {sa_lambda:.4f} mm")
        print(f"  SA-inf:    {sa_inf:.4f}")
        return sa_lambda, sa_inf, bin_centers, bin_means, \
               exponential_decay(bin_centers, *popt)
    except Exception as e:
        print(f"  Curve fitting failed: {e}")
        return np.nan, np.nan, bin_centers, bin_means, None


# =============================================================================
# PROCESS ONE SCAN
# =============================================================================

def process_scan(ts_path, fc_path, atlas_data, atlas_affine,
                 isocortex_labels, subcortex_labels,
                 scan_out, sub, ses, out_run, cohort=None):

    # Extract rectangular cross-block
    fc_cross, ts_combined, kept_iso, kept_sub, kept_all = filter_fc_cross(
        fc_path, ts_path, isocortex_labels, subcortex_labels
    )
    print(f"    Iso regions: {len(kept_iso)} | Sub regions: {len(kept_sub)}")
    print(f"    Cross-block shape: {fc_cross.shape}")

    # Centroids for iso and sub separately (needed for cross-distance matrix)
    iso_centroids = extract_centroids(kept_iso, atlas_data, atlas_affine)
    sub_centroids = extract_centroids(kept_sub, atlas_data, atlas_affine)

    # TA — computed across all iso+sub regions combined
    global_ta, ta_df = calculate_ta(ts_combined, kept_all)
    all_centroids = np.vstack([iso_centroids, sub_centroids])
    ta_df['x'] = all_centroids[:, 0]
    ta_df['y'] = all_centroids[:, 1]
    ta_df['z'] = all_centroids[:, 2]
    # Tag each region as iso or sub for reference
    ta_df['region_type'] = (
        ['isocortex'] * len(kept_iso) + ['subcortex'] * len(kept_sub)
    )
    ta_df.to_csv(os.path.join(scan_out, 'ta_results.csv'), index=False)

    # SA — computed on the cross-block with cross-distances
    sa_lambda, sa_inf, bin_x, bin_y, fit_y = calculate_sa_cross(
        fc_cross, iso_centroids, sub_centroids
    )
    if bin_x is not None and len(bin_x) > 0:
        pd.DataFrame({
            'Distance_mm':      bin_x,
            'Mean_Correlation': bin_y,
            'Fitted_Model':     fit_y if fit_y is not None
                                else [np.nan] * len(bin_x)
        }).to_csv(os.path.join(scan_out, 'sa_curve_fit.csv'), index=False)

    result = {
        'Subject':     sub,
        'Session':     ses,
        'Run':         out_run,
        'Global_TA':   global_ta,
        'SA_Lambda':   sa_lambda,
        'SA_Infinity': sa_inf
    }
    if cohort:
        result['Cohort'] = cohort
    print(f"    ✅  TA={global_ta:.4f}  "
          f"SA_lambda={sa_lambda:.4f}  SA_inf={sa_inf:.4f}")
    return result


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Batch TA/SA — Iso-Subcortex cross edges only"
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--parent_dir',
                      help='Cohort mode: parent folder containing cohort subdirectories')
    mode.add_argument('--timeseries_dir',
                      help='Direct mode: flat folder of timeseries CSVs')

    parser.add_argument('--fc_dir',
                        help='Direct mode: flat folder of FC matrix subfolders')
    parser.add_argument('--atlas',
                        help='Direct mode: path to atlas NIfTI (.nii.gz)')
    parser.add_argument('--hierarchy', required=True,
                        help='Path to ambaCCFv3_atlas_annotationsTable.csv')
    parser.add_argument('--output_dir', required=True,
                        help='Directory to save results')
    parser.add_argument('--ignore_direction', action='store_true',
                        help='Direct mode: ignore dir-AP/PA in filenames')
    parser.add_argument('--exclude_subjects', nargs='+', default=[])
    parser.add_argument('--exclude_sessions',  nargs='+', default=[])
    parser.add_argument('--exclude_runs',      nargs='+', default=[])
    parser.add_argument('--exclude_scans',     nargs='+', default=[])
    args = parser.parse_args()

    if args.timeseries_dir and (not args.fc_dir or not args.atlas):
        parser.error('--timeseries_dir requires --fc_dir and --atlas')

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading region label definitions from hierarchy...")
    isocortex_labels = get_isocortex_labels(args.hierarchy)
    subcortex_labels = get_subcortex_labels(args.hierarchy)
    print(f"  Isocortex labels: {len(isocortex_labels)}")
    print(f"  Subcortex labels: {len(subcortex_labels)}")

    summary_results = []

    # =========================================================================
    # DIRECT MODE
    # =========================================================================
    if args.timeseries_dir:
        print(f"\n{'='*60}")
        print("DIRECT MODE")
        print(f"{'='*60}")

        atlas_img    = nib.load(args.atlas)
        atlas_data   = atlas_img.get_fdata()
        atlas_affine = atlas_img.affine
        print(f"  Atlas: {os.path.basename(args.atlas)}, shape: {atlas_data.shape}")

        ts_files = sorted(glob.glob(
            os.path.join(args.timeseries_dir, '*.csv')))
        print(f"  Found {len(ts_files)} timeseries files")

        for ts_path in ts_files:
            fname = os.path.splitext(os.path.basename(ts_path))[0]
            sub, ses, direction, run = parse_bids(
                fname, has_direction=not args.ignore_direction)
            out_run = output_run_label(direction, run)
            scan_id = f"{sub}_{ses}_{out_run}"

            print(f"\n  Processing: {sub} | {ses} | {out_run}")

            if sub in args.exclude_subjects:
                print("    Skipping: subject excluded"); continue
            if ses in args.exclude_sessions:
                print("    Skipping: session excluded"); continue
            if out_run in args.exclude_runs:
                print("    Skipping: run excluded"); continue
            if scan_id in args.exclude_scans:
                print("    Skipping: scan excluded"); continue

            fc_path = find_fc_file(args.fc_dir, sub, ses, direction, run)
            if not fc_path:
                print("    Skipping: no matching FC matrix found"); continue

            print(f"    TS: {os.path.basename(ts_path)}")
            print(f"    FC: {os.path.basename(fc_path)}")

            scan_out = os.path.join(args.output_dir, scan_id)
            os.makedirs(scan_out, exist_ok=True)

            try:
                result = process_scan(
                    ts_path, fc_path, atlas_data, atlas_affine,
                    isocortex_labels, subcortex_labels,
                    scan_out, sub, ses, out_run
                )
                summary_results.append(result)
            except Exception as e:
                print(f"    ❌ Error: {e}")
                import traceback; traceback.print_exc()

    # =========================================================================
    # COHORT MODE
    # =========================================================================
    else:
        for cohort in COHORT_NAMES:
            cohort_dir = os.path.join(args.parent_dir, cohort)
            if not os.path.exists(cohort_dir):
                print(f"\nSkipping cohort '{cohort}' — directory not found")
                continue

            ts_root = os.path.join(cohort_dir, 'timeseries_csv')
            fc_root = os.path.join(cohort_dir, 'fc_matrix')

            if not os.path.exists(ts_root):
                print(f"\nSkipping '{cohort}' — timeseries_csv not found"); continue
            if not os.path.exists(fc_root):
                print(f"\nSkipping '{cohort}' — fc_matrix not found"); continue

            print(f"\n{'='*60}")
            print(f"COHORT: {cohort}")
            print(f"{'='*60}")

            atlas_files = glob.glob(os.path.join(cohort_dir, '*.nii.gz'))
            if len(atlas_files) == 0:
                print("  Skipping — no .nii.gz atlas found"); continue
            if len(atlas_files) > 1:
                print("  Skipping — multiple .nii.gz files found"); continue

            atlas_img    = nib.load(atlas_files[0])
            atlas_data   = atlas_img.get_fdata()
            atlas_affine = atlas_img.affine
            print(f"  Atlas: {os.path.basename(atlas_files[0])}, "
                  f"shape: {atlas_data.shape}")

            ts_files      = sorted(glob.glob(os.path.join(ts_root, '*.csv')))
            has_direction = cohort in COHORTS_WITH_DIRECTION
            print(f"  Found {len(ts_files)} timeseries files")

            for ts_path in ts_files:
                fname = os.path.splitext(os.path.basename(ts_path))[0]
                sub, ses, direction, run = parse_bids(
                    fname, has_direction=has_direction)
                out_run = output_run_label(direction, run, cohort)
                scan_id = f"{sub}_{ses}_{out_run}"

                print(f"\n  Processing: {sub} | {ses} | "
                      + (f"dir-{direction} | " if has_direction else "")
                      + f"run-{run} → {out_run}")

                if sub in args.exclude_subjects:
                    print("    Skipping: subject excluded"); continue
                if ses in args.exclude_sessions:
                    print("    Skipping: session excluded"); continue
                if f'run-{run}' in args.exclude_runs:
                    print("    Skipping: run excluded"); continue
                if scan_id in args.exclude_scans:
                    print("    Skipping: scan excluded"); continue

                fc_path = find_fc_file(fc_root, sub, ses, direction, run)
                if not fc_path:
                    print("    Skipping: no matching FC matrix found"); continue

                print(f"    TS: {os.path.basename(ts_path)}")
                print(f"    FC: {os.path.basename(fc_path)}")

                scan_out = os.path.join(args.output_dir, cohort, scan_id)
                os.makedirs(scan_out, exist_ok=True)

                try:
                    result = process_scan(
                        ts_path, fc_path, atlas_data, atlas_affine,
                        isocortex_labels, subcortex_labels,
                        scan_out, sub, ses, out_run, cohort
                    )
                    summary_results.append(result)
                except Exception as e:
                    print(f"    ❌ Error: {e}")
                    import traceback; traceback.print_exc()

    # =========================================================================
    # MASTER SUMMARY
    # =========================================================================
    if summary_results:
        print("\nWriting master summary...")
        summary_df = pd.DataFrame(summary_results)
        summary_df = sort_summary(summary_df)
        cols = [c for c in ['Cohort', 'Subject', 'Session', 'Run',
                             'Global_TA', 'SA_Lambda', 'SA_Infinity']
                if c in summary_df.columns]
        summary_df = summary_df[cols]
        out_path = os.path.join(
            args.output_dir,
            'all_subjects_isosubcortex_autocorrelation.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Done! Results saved to: {out_path}")
    else:
        print("\nNo scans were successfully processed.")


if __name__ == '__main__':
    main()
