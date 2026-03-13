"""
region_filter.py

Shared utility for filtering FC matrices and timeseries CSVs to a
specific set of brain region label values, based on the Allen CCF v3
ontology hierarchy.

Used by:
    batch_autocorrelation_isocortex.py
    batch_autocorrelation_subcortex.py
"""

import csv
import numpy as np
import pandas as pd


def get_isocortex_labels(hierarchy_csv):
    """
    Returns set of int label values whose structure_id_path contains '/315/'
    (315 = Isocortex node in Allen CCF v3). Includes both hemispheres.
    """
    labels = set()
    with open(hierarchy_csv, newline='') as f:
        for row in csv.DictReader(f):
            try:
                lid = int(row['labelID'])
            except (ValueError, KeyError):
                continue
            if '/315/' in row.get('structure_id_path', ''):
                labels.add(lid)
    return labels


def get_subcortex_labels(hierarchy_csv):
    """
    Returns set of int label values whose structure_id_path does NOT contain
    '/315/' — i.e. all non-isocortical grey matter regions.
    Excludes background (labelID <= 0).
    """
    all_labels = set()
    isocortex_labels = set()
    with open(hierarchy_csv, newline='') as f:
        for row in csv.DictReader(f):
            try:
                lid = int(row['labelID'])
            except (ValueError, KeyError):
                continue
            if lid <= 0:
                continue
            all_labels.add(lid)
            if '/315/' in row.get('structure_id_path', ''):
                isocortex_labels.add(lid)
    return all_labels - isocortex_labels


def filter_fc_cross(fc_path, ts_path, labels_rows, labels_cols):
    """
    Extracts a rectangular cross-block from the FC matrix:
        rows = regions in labels_rows (e.g. isocortex)
        cols = regions in labels_cols (e.g. subcortex)

    This gives the off-diagonal block of the full FC matrix representing
    connectivity between two distinct region sets (e.g. iso-subcortex edges).

    For TA, returns timeseries for ALL regions in both sets combined,
    since TA is computed per-region independently of the FC block.

    Returns:
        fc_cross      : np.ndarray (n_iso x n_sub) rectangular FC block
        ts_combined   : np.ndarray (n_timepoints x (n_iso + n_sub))
        kept_rows     : np.ndarray of int label values for rows (iso)
        kept_cols     : np.ndarray of int label values for cols (sub)
        kept_all      : np.ndarray of int label values for all regions (iso+sub)
    """
    # Load FC matrix
    fc_df = pd.read_csv(fc_path, header=0, index_col=0)
    if fc_df.isnull().values.any():
        fc_df = fc_df.fillna(0)
    fc_labels = fc_df.index.astype(float).astype(int).values

    # Load timeseries
    ts_raw = pd.read_csv(ts_path, header=None)
    ts_labels = ts_raw.iloc[0, :].values.astype(float).astype(int)
    ts_data_full = ts_raw.iloc[1:, :].values.astype(float)

    fc_set = set(fc_labels)
    ts_set = set(ts_labels)

    # Find labels present in both files for each region set
    common_rows = sorted(fc_set & ts_set & labels_rows)
    common_cols = sorted(fc_set & ts_set & labels_cols)

    if len(common_rows) < 1 or len(common_cols) < 1:
        raise ValueError(
            f"Insufficient cross-block regions: "
            f"{len(common_rows)} row regions, {len(common_cols)} col regions. "
            f"Need at least 1 of each."
        )

    # Preserve FC matrix ordering
    kept_rows = np.array([l for l in fc_labels if l in set(common_rows)])
    kept_cols = np.array([l for l in fc_labels if l in set(common_cols)])

    # Extract rectangular cross-block from FC matrix
    row_mask = np.isin(fc_labels, kept_rows)
    col_mask = np.isin(fc_labels, kept_cols)
    fc_cross  = fc_df.values[np.ix_(row_mask, col_mask)]

    # Build combined timeseries for all regions (rows + cols)
    kept_all = np.array([l for l in fc_labels
                         if l in set(common_rows) | set(common_cols)])
    label_to_ts_col = {int(ts_labels[i]): i for i in range(len(ts_labels))}
    ts_col_order    = [label_to_ts_col[l] for l in kept_all]
    ts_combined     = ts_data_full[:, ts_col_order]

    return fc_cross, ts_combined, kept_rows, kept_cols, kept_all


def filter_fc_and_ts(fc_path, ts_path, target_labels):
    """
    Loads and filters the FC matrix and timeseries to only the regions
    in target_labels that are actually present in both files.

    FC matrix CSV:
        - Row 0 / Col 0 are label intensity values (index + header)
        - Values are correlation coefficients

    Timeseries CSV:
        - Row 0 contains label intensity values (no pandas header)
        - Remaining rows are timepoints x regions

    Args:
        fc_path       : path to *fc_matrix.csv
        ts_path       : path to timeseries CSV
        target_labels : set of int label values to keep

    Returns:
        fc_matrix     : np.ndarray (n_regions x n_regions)
        ts_data       : np.ndarray (n_timepoints x n_regions)
        kept_labels   : np.ndarray of int label values retained
        
    Raises:
        ValueError if fewer than 2 matching regions are found
    """
    # --- Load FC matrix ---
    fc_df = pd.read_csv(fc_path, header=0, index_col=0)
    if fc_df.isnull().values.any():
        fc_df = fc_df.fillna(0)

    fc_labels = fc_df.index.astype(float).astype(int).values

    # --- Load timeseries ---
    ts_raw = pd.read_csv(ts_path, header=None)
    ts_labels = ts_raw.iloc[0, :].values.astype(float).astype(int)
    ts_data_full = ts_raw.iloc[1:, :].values.astype(float)

    # --- Find labels present in BOTH files AND in target set ---
    fc_set = set(fc_labels)
    ts_set = set(ts_labels)
    common = (fc_set & ts_set & target_labels)

    if len(common) < 2:
        raise ValueError(
            f"Only {len(common)} matching target regions found across FC matrix "
            f"and timeseries. Need at least 2 to run analysis."
        )

    # Preserve original ordering from FC matrix
    kept_labels = np.array([l for l in fc_labels if l in common])

    # Filter FC matrix
    keep_fc = np.isin(fc_labels, kept_labels)
    fc_filtered = fc_df.values[np.ix_(keep_fc, keep_fc)]

    # Filter timeseries (reorder to match FC label order)
    label_to_ts_col = {int(ts_labels[i]): i for i in range(len(ts_labels))}
    ts_col_order = [label_to_ts_col[l] for l in kept_labels]
    ts_filtered = ts_data_full[:, ts_col_order]

    return fc_filtered, ts_filtered, kept_labels
