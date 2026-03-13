"""
batch_autocorrelation_wholebrain.py

Batch spatial and temporal autocorrelation — WHOLE BRAIN (no region filtering).

Loops over all 5 cohort folders under a single parent directory:
    cohort_100, cohort_200_AP, cohort_200_PA, cohort_300_AP, cohort_300_PA

Within each cohort, expects:
    <cohort>/timeseries_csv/   — one CSV per scan, filename contains BIDS fields
    <cohort>/fc_matrix/        — one subfolder per scan containing *fc_matrix.csv
    <cohort>/*.nii.gz          — exactly one atlas NIfTI (auto-detected per cohort)

No region filtering is applied — all regions in the FC matrix and timeseries
are used as-is.

cohort_100: two real PA runs — run numbers kept as-is (run-1 and run-2).
cohort_200/300: PA acquisitions relabelled as run-2 in output.

Usage:
    python batch_autocorrelation_wholebrain.py \
        --parent_dir /path/to/cohorts_parent \
        --output_dir /path/to/output

    # Optional exclusions:
        --exclude_subjects sub-01 sub-02 \
        --exclude_sessions ses-06m \
        --exclude_runs run-1 \
        --exclude_scans sub-01_ses-06m_run-1
"""

import os
import argparse
import glob
import re
import numpy as np
import pandas as pd
import nibabel as nib

from calculate_autocorrelation import extract_centroids, calculate_ta, calculate_sa


# =============================================================================
# COHORT CONFIGURATION
# =============================================================================
# Cohorts are specified at runtime via --cohorts or --auto_cohorts.
# There is no hardcoded default list.


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


def output_run_label(direction, run, cohort, ignore_direction=False):
    """
    ignore_direction=True: always keep original run number as-is (for datasets
        with multiple real runs that don't need AP/PA remapping).
    cohort_100: two real PA runs — keep original run numbers as-is.
    cohort_200/300: PA is always a single acquisition → relabel as run-2.
    AP acquisitions always keep their original run number.
    """
    if ignore_direction:
        return f'run-{run}' if not str(run).startswith('run-') else run
    if direction == 'PA' and cohort != 'cohort_100':
        return 'run-2'
    return f'run-{run}' if not str(run).startswith('run-') else run


# =============================================================================
# FC MATRIX LOADER (no filtering)
# =============================================================================

def find_fc_file(fc_root, sub, ses, direction, run, ignore_direction=False):
    if not os.path.exists(fc_root):
        return None
    # Ensure run token is in "run-N" form for exact matching
    run_token = f'run-{run}' if not str(run).startswith('run-') else run
    for d in os.listdir(fc_root):
        dpath = os.path.join(fc_root, d)
        if not os.path.isdir(dpath):
            continue
        if sub not in d or ses not in d or run_token not in d:
            continue
        if not ignore_direction and direction != 'unknown' and f'dir-{direction}' not in d:
            continue
        csv_files = glob.glob(os.path.join(dpath, '*fc_matrix.csv'))
        if csv_files:
            return csv_files[0]
    return None


def load_fc_and_ts(fc_path, ts_path):
    """
    Loads FC matrix and timeseries without any region filtering.

    FC matrix: rows/cols are region labels (first row and first col are headers).
    Timeseries: first row = region labels, remaining rows = timepoints x regions.

    Returns:
        fc_matrix  : np.ndarray (n_regions x n_regions)
        ts_data    : np.ndarray (n_timepoints x n_regions)
        labels     : list of int region labels
    """
    # --- FC matrix ---
    fc_df = pd.read_csv(fc_path, index_col=0)
    fc_df.columns = fc_df.columns.astype(str)
    fc_df.index   = fc_df.index.astype(str)

    col_labels = [int(float(c)) for c in fc_df.columns]
    row_labels = [int(float(r)) for r in fc_df.index]

    # Use column labels as the canonical label set
    labels    = col_labels
    fc_matrix = fc_df.values.astype(float)

    # --- Timeseries ---
    ts_raw = pd.read_csv(ts_path, header=None)
    ts_labels = [int(float(v)) for v in ts_raw.iloc[0].tolist()]
    ts_data   = ts_raw.iloc[1:].values.astype(float)

    # Align timeseries columns to FC label order
    label_to_idx = {lbl: i for i, lbl in enumerate(ts_labels)}
    col_indices  = [label_to_idx[lbl] for lbl in labels if lbl in label_to_idx]
    labels       = [labels[i] for i in range(len(labels)) if labels[i] in label_to_idx]
    fc_matrix    = fc_matrix[np.ix_(
        [row_labels.index(l) for l in labels],
        [col_labels.index(l) for l in labels]
    )]
    ts_data = ts_data[:, col_indices]

    return fc_matrix, ts_data, labels


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
# MAIN
# =============================================================================

def process_dataset(ts_root, fc_root, atlas_path, output_dir, label,
                    ignore_direction, exclude_subjects, exclude_sessions,
                    exclude_runs, exclude_scans):
    print(f"\n{'='*60}")
    print(f"DATASET: {label}")
    print(f"  timeseries_csv : {ts_root}")
    print(f"  fc_matrix      : {fc_root}")
    print(f"  atlas          : {atlas_path}")
    print(f"{'='*60}")

    atlas_img    = nib.load(atlas_path)
    atlas_data   = atlas_img.get_fdata()
    atlas_affine = atlas_img.affine
    print(f"  Atlas shape: {atlas_data.shape}")

    ts_files      = sorted(glob.glob(os.path.join(ts_root, '*.csv')))
    has_direction = not ignore_direction
    print(f"  Found {len(ts_files)} timeseries files")

    summary_results = []

    for ts_path in ts_files:
        fname   = os.path.splitext(os.path.basename(ts_path))[0]
        sub, ses, direction, run = parse_bids(fname, has_direction=has_direction)
        out_run = output_run_label(direction, run, label,
                                   ignore_direction=ignore_direction)
        scan_id = f"{sub}_{ses}_{out_run}"

        print(f"\n  Processing: {sub} | {ses} | "
              + (f"dir-{direction} | " if has_direction else "")
              + f"run-{run} → {out_run}")

        if sub in exclude_subjects:
            print(f"    Skipping: subject excluded"); continue
        if ses in exclude_sessions:
            print(f"    Skipping: session excluded"); continue
        if f'run-{run}' in exclude_runs:
            print(f"    Skipping: run excluded"); continue
        if scan_id in exclude_scans:
            print(f"    Skipping: scan excluded"); continue

        fc_path = find_fc_file(fc_root, sub, ses, direction, run,
                                ignore_direction=ignore_direction)
        if not fc_path:
            print(f"    Skipping: no matching FC matrix found")
            continue

        print(f"    TS:  {os.path.basename(ts_path)}")
        print(f"    FC:  {os.path.basename(fc_path)}")

        scan_out = os.path.join(output_dir, label, scan_id)
        os.makedirs(scan_out, exist_ok=True)

        try:
            fc_matrix, ts_data, labels = load_fc_and_ts(fc_path, ts_path)
            print(f"    Regions loaded: {len(labels)}")

            centroids = extract_centroids(labels, atlas_data, atlas_affine)

            global_ta, ta_df = calculate_ta(ts_data, labels)
            ta_df['x'] = centroids[:, 0]
            ta_df['y'] = centroids[:, 1]
            ta_df['z'] = centroids[:, 2]
            ta_df.to_csv(os.path.join(scan_out, 'ta_results.csv'), index=False)

            sa_lambda, sa_inf, bin_x, bin_y, fit_y = calculate_sa(fc_matrix, centroids)
            if bin_x is not None and len(bin_x) > 0:
                pd.DataFrame({
                    'Distance_mm':      bin_x,
                    'Mean_Correlation': bin_y,
                    'Fitted_Model':     fit_y if fit_y is not None
                                        else [np.nan] * len(bin_x)
                }).to_csv(os.path.join(scan_out, 'sa_curve_fit.csv'), index=False)

            summary_results.append({
                'Cohort':      label,
                'Subject':     sub,
                'Session':     ses,
                'Run':         out_run,
                'Global_TA':   global_ta,
                'SA_Lambda':   sa_lambda,
                'SA_Infinity': sa_inf
            })
            print(f"    ✅  TA={global_ta:.4f}  "
                  f"SA_lambda={sa_lambda:.4f}  SA_inf={sa_inf:.4f}")

        except Exception as e:
            print(f"    ❌ Error: {e}")
            import traceback; traceback.print_exc()

    return summary_results


def main():
    parser = argparse.ArgumentParser(
        description='Batch TA/SA — Whole brain, no region filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:

  1) Direct (flat folder — timeseries_csv, fc_matrix, atlas all in one place):
       --timeseries_dir /path/to/timeseries_csv
       --fc_dir         /path/to/fc_matrix
       --atlas          /path/to/atlas.nii.gz
       --label          my_dataset

  2) Cohort (multi-cohort subfolders):
       --parent_dir /path/to/parent
       --cohorts cohort_A cohort_B   OR   --auto_cohorts

Example (direct mode):
    python batch_autocorrelation_wholebrain.py \\
        --timeseries_dir /path/iso+subcortex/timeseries_csv \\
        --fc_dir         /path/iso+subcortex/fc_matrix \\
        --atlas          /path/iso+subcortex/atlas.nii.gz \\
        --label          iso+subcortex \\
        --output_dir     /path/wholebrain_analysis \\
        --ignore_direction
        """
    )

    # Direct mode
    parser.add_argument('--timeseries_dir',
                        help='[Direct] Path to timeseries_csv folder')
    parser.add_argument('--fc_dir',
                        help='[Direct] Path to fc_matrix folder')
    parser.add_argument('--atlas',
                        help='[Direct] Path to atlas .nii.gz')
    parser.add_argument('--label', default='dataset',
                        help='[Direct] Label used as subfolder name in output')

    # Cohort mode
    parser.add_argument('--parent_dir',
                        help='[Cohort] Parent folder containing cohort subdirs')
    parser.add_argument('--cohorts', nargs='+', default=None,
                        help='[Cohort] Cohort subdirectory names')
    parser.add_argument('--auto_cohorts', action='store_true',
                        help='[Cohort] Auto-detect cohort subdirs in parent_dir')

    # Shared
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--ignore_direction', action='store_true',
                        help='Ignore AP/PA direction — keep original run numbers as-is')
    parser.add_argument('--exclude_subjects', nargs='+', default=[])
    parser.add_argument('--exclude_sessions',  nargs='+', default=[])
    parser.add_argument('--exclude_runs',      nargs='+', default=[])
    parser.add_argument('--exclude_scans',     nargs='+', default=[])
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    direct_mode = bool(args.timeseries_dir or args.fc_dir or args.atlas)
    cohort_mode = bool(args.parent_dir or args.cohorts or args.auto_cohorts)

    if direct_mode and cohort_mode:
        parser.error("Cannot mix --timeseries_dir/--fc_dir/--atlas "
                     "with --parent_dir/--cohorts/--auto_cohorts.")

    summary_results = []

    # --- Direct mode ---
    if direct_mode:
        for flag, val in [('--timeseries_dir', args.timeseries_dir),
                          ('--fc_dir',         args.fc_dir),
                          ('--atlas',          args.atlas)]:
            if not val:
                parser.error(f"Direct mode requires {flag}")
        summary_results = process_dataset(
            ts_root=args.timeseries_dir, fc_root=args.fc_dir,
            atlas_path=args.atlas, output_dir=args.output_dir,
            label=args.label, ignore_direction=args.ignore_direction,
            exclude_subjects=args.exclude_subjects,
            exclude_sessions=args.exclude_sessions,
            exclude_runs=args.exclude_runs, exclude_scans=args.exclude_scans,
        )

    # --- Cohort mode ---
    else:
        if not args.parent_dir:
            parser.error("Cohort mode requires --parent_dir")
        if args.auto_cohorts:
            cohort_names = sorted([
                d for d in os.listdir(args.parent_dir)
                if os.path.isdir(os.path.join(args.parent_dir, d))
                and os.path.exists(os.path.join(args.parent_dir, d, 'timeseries_csv'))
            ])
            print(f"Auto-detected {len(cohort_names)} cohort(s): {cohort_names}")
        elif args.cohorts:
            cohort_names = args.cohorts
        else:
            parser.error("Cohort mode requires --cohorts or --auto_cohorts")

        for cohort in cohort_names:
            cohort_dir  = cohort if os.path.isabs(cohort) \
                          else os.path.join(args.parent_dir, cohort)
            cohort_name = os.path.basename(cohort.rstrip('/'))

            if not os.path.exists(cohort_dir):
                print(f"\nSkipping '{cohort_name}' — not found: {cohort_dir}"); continue

            ts_root    = os.path.join(cohort_dir, 'timeseries_csv')
            fc_root    = os.path.join(cohort_dir, 'fc_matrix')
            atlas_list = glob.glob(os.path.join(cohort_dir, '*.nii.gz'))

            if not os.path.exists(ts_root):
                print(f"\nSkipping '{cohort_name}' — timeseries_csv not found"); continue
            if not os.path.exists(fc_root):
                print(f"\nSkipping '{cohort_name}' — fc_matrix not found"); continue
            if len(atlas_list) == 0:
                print(f"\nSkipping '{cohort_name}' — no .nii.gz atlas found"); continue
            if len(atlas_list) > 1:
                print(f"\nSkipping '{cohort_name}' — multiple .nii.gz found:"); continue

            summary_results.extend(process_dataset(
                ts_root=ts_root, fc_root=fc_root, atlas_path=atlas_list[0],
                output_dir=args.output_dir, label=cohort_name,
                ignore_direction=args.ignore_direction,
                exclude_subjects=args.exclude_subjects,
                exclude_sessions=args.exclude_sessions,
                exclude_runs=args.exclude_runs, exclude_scans=args.exclude_scans,
            ))

    # --- Master summary ---
    if summary_results:
        print("\nWriting master summary...")
        summary_df = pd.DataFrame(summary_results)
        summary_df = sort_summary(summary_df)
        summary_df = summary_df[['Cohort', 'Subject', 'Session', 'Run',
                                  'Global_TA', 'SA_Lambda', 'SA_Infinity']]
        out_path = os.path.join(
            args.output_dir, 'all_subjects_wholebrain_autocorrelation.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Done! Results saved to: {out_path}")
    else:
        print("\nNo scans were successfully processed.")


if __name__ == '__main__':
    main()
