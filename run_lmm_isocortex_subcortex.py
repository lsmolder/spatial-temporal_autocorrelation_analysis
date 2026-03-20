#!/usr/bin/env python3
"""
run_lmm_isocortex_subcortex.py

Linear mixed model analysis of SA and TA metrics for one or more region sets
(isocortex, subcortex, wholebrain), followed by parcel-level TA driver analysis
and distance-dependent FC plots.

Key features
------------
- --filter_genotypes   : restrict analysis to specified genotype(s)
- --outlier_scope      : 'within' (per-subject MAD across runs) or
                         'between' (cross-subject MAD on run-pair delta)
                         Outlier removal is per-metric and independent.
- Runs are averaged within Subject x Session before the LMM.
- Distance x FC plot produced per genotype group for each analysis.

Usage
-----
    python run_lmm_isocortex_subcortex.py \
        --analyses wholebrain \
        --csv  /path/wholebrain.csv \
        --batch /path/wb_batch \
        --participants /path/participants.tsv \
        --model_vars Genotype Session \
        --outlier_scope within \
        --filter_genotypes KI2APOE4 KI2APOE3 \
        --n_dist_bins 20 \
        --output_dir /path/output
"""

import argparse
import os
import glob
import re
import warnings
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ============================================================
# Constants
# ============================================================

METRICS          = ["Global_TA", "SA_Lambda", "SA_Infinity"]
MAD_MULTIPLIER   = 3.0

REQUIRED_COLUMNS = ['Subject', 'Session', 'Run', 'Treatment', 'Genotype', 'sex']

PARTICIPANT_COLS = {
    'Treatment': {"drug": "Drug", "Drug": "Drug",
                  "vehicle": "Vehicle", "Vehicle": "Vehicle"},
    'Genotype':  {}
}

REFERENCE_LEVELS = {
    'Treatment': 'Vehicle',
    'Genotype':  None,
    'Session':   None,
}


# ============================================================
# Helpers
# ============================================================

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def mad(arr: np.ndarray) -> float:
    """Median absolute deviation (no normalisation factor)."""
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return np.nan
    return float(np.median(np.abs(arr - np.median(arr))))


def standardize_col(series: pd.Series, col: str) -> pd.Series:
    mapping = PARTICIPANT_COLS.get(col, {})
    if mapping:
        return series.astype(str).map(mapping).fillna(series.astype(str))
    return series.astype(str)


def normalize_subject(s: str) -> str:
    if s.lower().startswith('sub-'):
        return 'sub-' + s[4:]
    return s


def ses_sort_key(s):
    d = re.findall(r'\d+', str(s))
    return int(d[0]) if d else float('inf')


def run_sort_key(r):
    d = re.findall(r'\d+', str(r))
    return int(d[0]) if d else float('inf')


def mean_se(x):
    arr = pd.to_numeric(x, errors='coerce').dropna().to_numpy()
    n = len(arr)
    if n == 0:
        return np.nan, np.nan, 0
    return arr.mean(), (arr.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0), n


# ============================================================
# Model formula builder
# ============================================================

def build_formula(metric: str, model_vars: list, ref_levels: dict) -> tuple:
    short_map     = {f'v{i}': v for i, v in enumerate(model_vars)}
    formula_parts = []
    for short, var in short_map.items():
        ref = ref_levels.get(var)
        if ref:
            formula_parts.append(f"C({short}, Treatment(reference='{ref}'))")
        else:
            formula_parts.append(f"C({short})")
    formula    = f"{metric} ~ " + " * ".join(formula_parts)
    rename_map = {}
    for short, var in short_map.items():
        ref = ref_levels.get(var)
        old = (f"C({short}, Treatment(reference='{ref}'))"
               if ref else f"C({short})")
        new = (f"C({var}, Treatment(reference='{ref}'))"
               if ref else f"C({var})")
        rename_map[old] = new
    return formula, short_map, rename_map


def rename_terms(terms: list, rename_map: dict) -> list:
    result = []
    for term in terms:
        for old, new in rename_map.items():
            term = term.replace(old, new)
        result.append(term)
    return result


# ============================================================
# Data loading
# ============================================================

def load_and_merge(metrics_csv: str, participants_tsv: str,
                   model_vars: list) -> pd.DataFrame:
    df = pd.read_csv(metrics_csv)
    df['Subject'] = df['Subject'].astype(str).str.strip().apply(normalize_subject)
    df['Session'] = df['Session'].astype(str)
    df['Run']     = df['Run'].astype(str)

    participants = pd.read_csv(participants_tsv, sep='\t')
    participants.columns = participants.columns.str.strip()  # remove leading/trailing spaces from column names
    # strip whitespace from all string columns
    for col in participants.select_dtypes(include='object').columns:
        participants[col] = participants[col].astype(str).str.strip()
    col_map = {c.lower(): c for c in participants.columns}

    if 'subject' in col_map and col_map['subject'] != 'Subject':
        participants = participants.rename(columns={col_map['subject']: 'Subject'})
    elif 'subject' not in col_map:
        raise ValueError(f"No 'Subject' column in participants TSV. "
                         f"Columns: {list(participants.columns)}")

    participants['Subject'] = (participants['Subject'].astype(str)
                                .str.strip().apply(normalize_subject))

    rename = {}
    for canonical in PARTICIPANT_COLS:
        lower = canonical.lower()
        if lower in col_map and col_map[lower] != canonical:
            rename[col_map[lower]] = canonical
    # Handle sex column capitalisation separately (stored as 'sex' not in PARTICIPANT_COLS)
    if 'sex' in col_map and col_map['sex'] != 'sex':
        rename[col_map['sex']] = 'sex'
    if rename:
        participants = participants.rename(columns=rename)

    for col in PARTICIPANT_COLS:
        if col in participants.columns:
            participants[col] = standardize_col(participants[col], col)

    MERGE_KEYS = ['Subject', 'Session', 'Run']
    for tsv_col, canonical in [('session', 'Session'), ('run', 'Run')]:
        if tsv_col in col_map and col_map[tsv_col] != canonical:
            participants = participants.rename(columns={col_map[tsv_col]: canonical})
        elif canonical in participants.columns and tsv_col in participants.columns:
            participants = participants.drop(columns=[tsv_col])

    if 'Session' in participants.columns:
        participants['Session'] = participants['Session'].astype(str).str.strip()
    if 'Run' in participants.columns:
        participants['Run'] = participants['Run'].astype(str).str.strip()

    available_keys = [k for k in MERGE_KEYS if k in participants.columns]
    print(f"  Merging on: {available_keys}")

    tsv_non_key = [c for c in participants.columns if c not in available_keys]
    conflict    = [c for c in tsv_non_key if c in df.columns]
    if conflict:
        participants = participants.drop(columns=conflict)
        print(f"  Dropped duplicate TSV columns: {conflict}")

    df = df.merge(participants, on=available_keys, how='left')

    for col in PARTICIPANT_COLS:
        if col not in df.columns:
            print(f"  Note: '{col}' not in TSV — column will be blank")
        else:
            n = df[col].isna().sum()
            if n > 0:
                print(f"  Note: {n} rows have no '{col}'")

    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            df[col] = ''
        df[col] = df[col].fillna('')

    metric_cols = [c for c in METRICS if c in df.columns]
    cohort_col  = ['Cohort'] if 'Cohort' in df.columns else []
    fixed_cols  = ['Subject', 'Session', 'Run'] + cohort_col + metric_cols + \
                  ['Treatment', 'Genotype', 'sex']
    extra_cols  = [c for c in df.columns if c not in fixed_cols]
    df = df[[c for c in fixed_cols if c in df.columns] + extra_cols]
    return df


# ============================================================
# Genotype filter
# ============================================================

def apply_genotype_filter(df: pd.DataFrame,
                           filter_genotypes: list) -> pd.DataFrame:
    if not filter_genotypes:
        return df
    before = len(df)
    df = df[df['Genotype'].astype(str).isin(filter_genotypes)].copy()
    print(f"  Genotype filter {filter_genotypes}: "
          f"{before} → {len(df)} rows "
          f"({df['Subject'].nunique()} subjects)")
    return df


# ============================================================
# Prepare data
# ============================================================

def prepare_data(df: pd.DataFrame, model_vars: list,
                 ref_levels: dict, outdir: Path) -> tuple:
    dat = df.copy()

    for var in model_vars:
        if var == 'Session':
            sessions = sorted(dat['Session'].dropna().unique(), key=ses_sort_key)
            if ref_levels.get('Session') is None:
                ref_levels['Session'] = sessions[0] if sessions else None
            dat['Session'] = pd.Categorical(dat['Session'],
                                            categories=sessions, ordered=True)
        elif var == 'Run':
            runs = sorted(dat['Run'].dropna().unique(), key=run_sort_key)
            dat['Run'] = pd.Categorical(dat['Run'], categories=runs, ordered=True)
        else:
            vals = sorted(dat[var].dropna().astype(str).unique())
            if ref_levels.get(var) is None:
                ref_levels[var] = vals[0] if vals else None
                print(f"  Reference level for '{var}': '{ref_levels[var]}'")
            dat[var] = pd.Categorical(dat[var], categories=vals, ordered=False)

    participant_vars = [v for v in model_vars if v not in ('Session', 'Run')]
    mask = pd.Series(True, index=dat.index)
    for var in participant_vars:
        mask = mask & (~dat[var].astype(str).isin(['NA', 'na', '']))
    n_excluded = (~mask).sum()
    if n_excluded > 0:
        print(f"  Excluding {n_excluded} rows where participant variables are blank")
    dat = dat[mask].copy()

    lines = [f"Total subjects: {dat['Subject'].nunique()}"]
    for var in model_vars:
        if var in dat.columns:
            lines.append(f"{var}: {sorted(dat[var].dropna().astype(str).unique())}")
    (outdir / 'subject_summary.txt').write_text('\n'.join(lines))
    return dat, ref_levels


# ============================================================
# Outlier removal — per metric, independent
# ============================================================

def remove_outliers_within(df: pd.DataFrame) -> tuple:
    """
    Within-subject MAD outlier removal.
    For each Subject x Session, compute median and MAD across their runs
    per metric. Flag any run outside median ± 3*MAD. Remove flagged runs.
    Handled independently per metric — a run flagged for Global_TA is only
    removed from Global_TA, not from SA_Lambda or SA_Infinity.
    """
    records = []
    df_clean = df.copy()

    for metric in METRICS:
        for (sub, ses), g in df.groupby(['Subject', 'Session'], observed=True):
            vals = g[metric].dropna()
            if len(vals) < 2:
                continue
            med   = float(np.median(vals))
            m     = mad(vals.to_numpy())
            if np.isnan(m) or m == 0:
                continue
            thresh = MAD_MULTIPLIER * m
            outlier_idx = g.index[
                (g[metric] - med).abs() > thresh
            ]
            for idx in outlier_idx:
                records.append({
                    'Subject':  sub,
                    'Session':  ses,
                    'Run':      df.loc[idx, 'Run'],
                    'Metric':   metric,
                    'Value':    df.loc[idx, metric],
                    'Median':   med,
                    'MAD':      m,
                    'Threshold':thresh,
                    'Scope':    'within'
                })
                df_clean.loc[idx, metric] = np.nan

    removed_df = pd.DataFrame(records)
    n = len(removed_df)
    print(f"  Within-subject outlier removal: {n} run-metric values set to NaN")
    if n > 0:
        for metric in METRICS:
            sub_df = removed_df[removed_df['Metric'] == metric]
            if not sub_df.empty:
                print(f"    {metric}: {len(sub_df)} removed")
    return df_clean, removed_df


def remove_outliers_between(df: pd.DataFrame) -> tuple:
    """
    Between-subject MAD outlier removal.
    For each Subject x Session, compute the range (max - min) across runs
    per metric. Collect all ranges across subjects. Flag any Subject x Session
    where range > median(ranges) + 3*MAD(ranges). Remove the single run
    furthest from that subject's mean for that metric.
    Handled independently per metric.
    """
    records   = []
    df_clean  = df.copy()

    for metric in METRICS:
        # Step 1: compute per-subject range
        ranges = {}
        for (sub, ses), g in df.groupby(['Subject', 'Session'], observed=True):
            vals = g[metric].dropna()
            if len(vals) < 2:
                continue
            ranges[(sub, ses)] = {
                'range': float(vals.max() - vals.min()),
                'g':     g
            }

        if not ranges:
            continue

        range_vals = np.array([v['range'] for v in ranges.values()])
        med_range  = float(np.median(range_vals))
        mad_range  = mad(range_vals)

        if np.isnan(mad_range) or mad_range == 0:
            continue

        thresh = med_range + MAD_MULTIPLIER * mad_range

        # Step 2: flag and remove worst run
        for (sub, ses), info in ranges.items():
            if info['range'] > thresh:
                g        = info['g']
                sub_mean = float(g[metric].mean())
                worst_idx = (g[metric] - sub_mean).abs().idxmax()
                records.append({
                    'Subject':       sub,
                    'Session':       ses,
                    'Run':           df.loc[worst_idx, 'Run'],
                    'Metric':        metric,
                    'Value':         df.loc[worst_idx, metric],
                    'Range':         info['range'],
                    'Median_Range':  med_range,
                    'MAD_Range':     mad_range,
                    'Threshold':     thresh,
                    'Scope':         'between'
                })
                df_clean.loc[worst_idx, metric] = np.nan

    removed_df = pd.DataFrame(records)
    n = len(removed_df)
    print(f"  Between-subject outlier removal: {n} run-metric values set to NaN")
    if n > 0:
        for metric in METRICS:
            sub_df = removed_df[removed_df['Metric'] == metric]
            if not sub_df.empty:
                print(f"    {metric}: {len(sub_df)} removed")
    return df_clean, removed_df


# ============================================================
# Average runs within Subject x Session
# ============================================================

def average_runs(df: pd.DataFrame, model_vars: list) -> pd.DataFrame:
    """
    Average all metric values within Subject x Session.
    This is the data that goes into the LMM.
    Non-metric columns (participant vars, etc.) are taken from the first row
    of each group since they are constant within subject x session.
    """
    group_cols = (['Subject', 'Session'] +
                  [v for v in model_vars
                   if v not in ('Session', 'Run') and v in df.columns])

    # Average metrics
    metric_avg = (df.groupby(group_cols, observed=True)[METRICS]
                  .mean().reset_index())

    # Grab non-metric, non-grouping extra cols from first row
    extra_cols = [c for c in df.columns
                  if c not in group_cols + METRICS + ['Run', 'Cohort']]
    if extra_cols:
        first_rows = (df.groupby(group_cols, observed=True)[extra_cols]
                      .first().reset_index())
        metric_avg = metric_avg.merge(first_rows, on=group_cols, how='left')

    if 'Cohort' in df.columns:
        cohort_first = (df.groupby(group_cols, observed=True)['Cohort']
                        .first().reset_index())
        metric_avg = metric_avg.merge(cohort_first, on=group_cols, how='left')

    n_runs_kept = df.groupby(group_cols, observed=True)['Run'].count().reset_index(
        name='n_runs_averaged')
    metric_avg = metric_avg.merge(n_runs_kept, on=group_cols, how='left')

    return metric_avg


# ============================================================
# Mixed model
# ============================================================

def fit_mixed_model(df: pd.DataFrame, metric: str, model_vars: list,
                    ref_levels: dict) -> pd.DataFrame:
    cols_needed = ['Subject'] + model_vars + [metric]
    sub = df[[c for c in cols_needed if c in df.columns]].copy()
    for var in [v for v in model_vars if v not in ('Session', 'Run')]:
        sub = sub[~sub[var].astype(str).isin(['NA', 'na', ''])]
    sub = sub.dropna(subset=[metric])

    empty = pd.DataFrame([{
        'metric': metric, 'term': np.nan,
        'estimate': np.nan, 'std_error': np.nan,
        'z_value': np.nan, 'p_value': np.nan,
        'n_obs': len(sub), 'n_subjects': sub['Subject'].nunique()
    }])

    if sub.empty or sub['Subject'].nunique() < 2:
        return empty

    formula, short_map, rename_map = build_formula(metric, model_vars, ref_levels)
    sub_model = sub.rename(columns={v: k for k, v in short_map.items()})

    try:
        result = smf.mixedlm(formula, data=sub_model,
                             groups=sub_model['Subject']).fit(
            reml=True, method='powell', maxiter=500, disp=False)
        terms = rename_terms(list(result.params.index), rename_map)
        out = pd.DataFrame({
            'term':      terms,
            'estimate':  result.params.values,
            'std_error': result.bse.values,
            'z_value':   result.tvalues.values,
            'p_value':   result.pvalues.values,
        })
        out['metric']     = metric
        out['n_obs']      = len(sub)
        out['n_subjects'] = sub['Subject'].nunique()
        return out[['metric', 'term', 'estimate', 'std_error',
                    'z_value', 'p_value', 'n_obs', 'n_subjects']]
    except Exception as e:
        return pd.DataFrame([{
            'metric': metric, 'term': f'MODEL_FAILED: {e}',
            'estimate': np.nan, 'std_error': np.nan,
            'z_value': np.nan, 'p_value': np.nan,
            'n_obs': len(sub), 'n_subjects': sub['Subject'].nunique()
        }])


def fit_all_models(df, model_vars, ref_levels, dataset_label):
    results = []
    for metric in METRICS:
        res = fit_mixed_model(df, metric, model_vars, ref_levels)
        res['dataset'] = dataset_label
        results.append(res)
    return pd.concat(results, ignore_index=True)


# ============================================================
# Plotting — LMM trajectory
# ============================================================

def plot_metric(avg_df, metric, model_vars, title_suffix, outbase):
    if avg_df.empty or metric not in avg_df.columns:
        return
    x_var    = 'Session' if 'Session' in model_vars else model_vars[0]
    line_var = next((v for v in model_vars if v != x_var), None)
    x_vals   = sorted(avg_df[x_var].dropna().astype(str).unique(),
                      key=ses_sort_key if x_var == 'Session' else str)
    x_to_int  = {v: i for i, v in enumerate(x_vals)}
    line_vals = (sorted(avg_df[line_var].dropna().astype(str).unique())
                 if line_var else [None])
    markers    = ['o', 's', '^', 'D', 'v', 'P']
    linestyles = ['-', '--', ':', '-.']

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for i, lv in enumerate(line_vals):
        sub = (avg_df[avg_df[line_var].astype(str) == lv]
               if line_var else avg_df)
        rows = []
        for xv in x_vals:
            m, se, _ = mean_se(sub[sub[x_var].astype(str) == xv][metric])
            rows.append({'x': x_to_int[xv], 'mean': m, 'se': se})
        s = pd.DataFrame(rows).dropna(subset=['mean'])
        ax.plot(s['x'], s['mean'],
                marker=markers[i % len(markers)],
                linestyle=linestyles[i % len(linestyles)],
                linewidth=1.5, markersize=6,
                label=str(lv) if line_var else metric)
        ax.errorbar(s['x'], s['mean'], yerr=s['se'],
                    fmt='none', capsize=3, linewidth=1.2)

    ax.set_xticks(range(len(x_vals)))
    ax.set_xticklabels(x_vals)
    ax.set_xlabel(x_var)
    ax.set_ylabel(metric)
    ax.set_title(f'{metric}{title_suffix}', fontweight='bold')
    if line_var:
        ax.legend(title=line_var, frameon=False)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(outbase.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)


def plot_raw_runlevel(dat, metric, model_vars, outbase):
    if 'OutlierRemoved' not in dat.columns:
        dat = dat.copy()
        dat['OutlierRemoved'] = 0

    x_var    = 'Session' if 'Session' in model_vars else model_vars[0]
    x_vals   = sorted(dat[x_var].dropna().astype(str).unique(),
                      key=ses_sort_key if x_var == 'Session' else str)
    x_to_int = {v: i for i, v in enumerate(x_vals)}
    line_var  = next((v for v in model_vars if v != x_var), None)
    line_vals = (sorted(dat[line_var].dropna().astype(str).unique())
                 if line_var else [None])
    markers   = ['o', 's', '^', 'D']

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for i, lv in enumerate(line_vals):
        sub    = (dat[dat[line_var].astype(str) == lv].copy()
                  if line_var else dat.copy())
        x      = sub[x_var].astype(str).map(x_to_int).astype(float).to_numpy()
        jitter = np.random.uniform(-0.06, 0.06, size=len(sub))
        y      = sub[metric].to_numpy()
        label  = str(lv) if line_var else metric
        kept   = ~np.isnan(y)
        ax.scatter(x[kept] + jitter[kept], y[kept],
                   marker=markers[i % len(markers)], alpha=0.7,
                   label=f'{label} (kept)')
        removed = np.isnan(y)
        if removed.any():
            ax.scatter(x[removed] + jitter[removed],
                       np.zeros(removed.sum()),
                       marker='x', alpha=0.5, color='red',
                       label=f'{label} (removed)')
        means   = sub.groupby(x_var, observed=True)[metric].mean()
        mean_xs = [x_to_int[str(v)] for v in means.index if str(v) in x_to_int]
        ax.plot(mean_xs, means.values,
                linestyle='-' if i == 0 else '--', linewidth=1.3)

    ax.set_xticks(range(len(x_vals)))
    ax.set_xticklabels(x_vals)
    ax.set_xlabel(x_var)
    ax.set_ylabel(metric)
    ax.set_title(f'{metric} raw run-level data', fontweight='bold')
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(outbase.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)


# ============================================================
# Distance x FC plot
# ============================================================

def _exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c


def collect_fc_matrices(batch_dir: str, summary_df: pd.DataFrame) -> dict:
    """
    Walk batch_dir for sa_curve_fit.csv files.
    Returns dict: {(Subject, Session, Run): (bin_centers, bin_means)}
    Also returns a subject->genotype map from summary_df.
    """
    curve_files = glob.glob(
        os.path.join(batch_dir, '**', 'sa_curve_fit.csv'), recursive=True)
    print(f"  Found {len(curve_files)} sa_curve_fit.csv files for distance-FC plot")

    sub_geno = {}
    if 'Genotype' in summary_df.columns:
        for _, row in summary_df.iterrows():
            sub_geno[str(row['Subject']).lower()] = str(row['Genotype'])

    data = {}
    for path in curve_files:
        folder = os.path.basename(os.path.dirname(path))
        sub_m  = re.search(r'(sub-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        ses_m  = re.search(r'(ses-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        run_m  = re.search(r'(run-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        if not (sub_m and ses_m and run_m):
            continue
        sub = sub_m.group(1)
        ses = ses_m.group(1)
        run = run_m.group(1)
        geno = sub_geno.get(sub.lower(), 'Unknown')
        try:
            curve_df = pd.read_csv(path)
            if 'Distance_mm' not in curve_df.columns or \
               'Mean_Correlation' not in curve_df.columns:
                continue
            data[(sub, ses, run, geno)] = (
                curve_df['Distance_mm'].to_numpy(),
                curve_df['Mean_Correlation'].to_numpy()
            )
        except Exception:
            continue
    return data


def plot_distance_fc_by_group(batch_dir: str, summary_df: pd.DataFrame,
                               outdir: Path, n_bins: int,
                               filter_genotypes: list) -> None:
    """
    Collect sa_curve_fit.csv files from batch_dir, bin distances,
    plot mean FC ± SEM per bin as points with exponential decay fit,
    one line per genotype group.
    """
    curve_data = collect_fc_matrices(batch_dir, summary_df)
    if not curve_data:
        print("  No sa_curve_fit.csv files found — skipping distance-FC plot")
        return

    # Group by genotype
    genotype_curves = {}
    for (sub, ses, run, geno), (dists, fcs) in curve_data.items():
        if filter_genotypes and geno not in filter_genotypes:
            continue
        if geno not in genotype_curves:
            genotype_curves[geno] = {'dists': [], 'fcs': []}
        genotype_curves[geno]['dists'].append(dists)
        genotype_curves[geno]['fcs'].append(fcs)

    if not genotype_curves:
        print("  No genotype groups found for distance-FC plot")
        return

    # Find global distance range
    all_dists = np.concatenate([
        d for v in genotype_curves.values() for d in v['dists']
    ])
    d_min, d_max = all_dists.min(), all_dists.max()
    bin_edges   = np.linspace(d_min, d_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    cmap   = plt.get_cmap('tab10')
    colors = [cmap(i) for i in range(len(genotype_curves))]

    fig, ax = plt.subplots(figsize=(8, 5))

    for idx, (geno, curves) in enumerate(sorted(genotype_curves.items())):
        color = colors[idx % len(colors)]

        # Pool all scan curves into per-bin arrays
        bin_fc_lists = [[] for _ in range(n_bins)]
        for dists, fcs in zip(curves['dists'], curves['fcs']):
            bin_idx = np.digitize(dists, bin_edges) - 1
            bin_idx = np.clip(bin_idx, 0, n_bins - 1)
            for b in range(n_bins):
                vals = fcs[bin_idx == b]
                if len(vals) > 0:
                    bin_fc_lists[b].append(float(np.mean(vals)))

        # Compute mean ± SEM per bin
        bx, by, bse = [], [], []
        for b, vals in enumerate(bin_fc_lists):
            if len(vals) == 0:
                continue
            arr = np.array(vals)
            bx.append(bin_centers[b])
            by.append(arr.mean())
            bse.append(arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0)

        bx  = np.array(bx)
        by  = np.array(by)
        bse = np.array(bse)

        # Plot binned points ± SEM
        ax.errorbar(bx, by, yerr=bse, fmt='o', color=color,
                    markersize=5, capsize=3, linewidth=1,
                    label=f'{geno} (binned mean ± SEM)')

        # Exponential decay fit
        if len(bx) >= 3:
            try:
                p0   = [by.max() - by.min(), 0.01, by.min()]
                popt, _ = curve_fit(_exp_decay, bx, by, p0=p0, maxfev=10000)
                x_fit = np.linspace(bx.min(), bx.max(), 200)
                ax.plot(x_fit, _exp_decay(x_fit, *popt),
                        color=color, linewidth=2, linestyle='--',
                        label=f'{geno} (exp fit)')
            except RuntimeError:
                pass

    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.set_xlabel('Inter-region Distance (mm)', fontsize=12)
    ax.set_ylabel('Mean Functional Connectivity (r)', fontsize=12)
    ax.set_title('Distance-Dependent FC by Genotype', fontweight='bold')
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / 'distance_fc_by_genotype.png', dpi=300,
                bbox_inches='tight')
    fig.savefig(outdir / 'distance_fc_by_genotype.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  Distance-FC plot saved → distance_fc_by_genotype.png")


# ============================================================
# Parcel-level TA driver analysis
# ============================================================

def collect_parcel_ta(batch_dir, summary_df, model_vars):
    print("  Collecting parcel-level TA from ta_results.csv files...")
    var_cols = [v for v in model_vars
                if v != 'Session' and v in summary_df.columns]

    # Build lookup on Subject + Session only — Run has been averaged away
    trt_map = {}
    for _, row in summary_df.iterrows():
        key = (str(row['Subject']).strip().lower(),
               str(row['Session']).strip())
        trt_map[key] = {v: row[v] for v in var_cols}

    ta_files = glob.glob(os.path.join(batch_dir, '**', 'ta_results.csv'),
                         recursive=True)
    print(f"  Found {len(ta_files)} ta_results.csv files")

    records = []
    for ta_path in ta_files:
        folder = os.path.basename(os.path.dirname(ta_path))
        sub_m  = re.search(r'(sub-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        ses_m  = re.search(r'(ses-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        run_m  = re.search(r'(run-[a-zA-Z0-9]+)', folder, re.IGNORECASE)
        if not (sub_m and ses_m and run_m):
            continue
        key      = (sub_m.group(1).lower(), ses_m.group(1))
        var_vals = trt_map.get(key)
        if var_vals is None:
            continue
        try:
            ta_df = pd.read_csv(ta_path)
        except Exception:
            continue
        if 'Region_Label' not in ta_df.columns or 'TA_Lag1' not in ta_df.columns:
            continue
        for _, r in ta_df.iterrows():
            rec = {'Subject': sub_m.group(1), 'Session': ses_m.group(1),
                   'Run': run_m.group(1),
                   'Region_Label': int(r['Region_Label']),
                   'TA_Lag1': float(r['TA_Lag1'])}
            rec.update(var_vals)
            records.append(rec)

    if not records:
        return pd.DataFrame()

    long_df = pd.DataFrame(records)

    # Average TA across runs within Subject x Session x Region
    group_cols = ['Subject', 'Session', 'Region_Label'] + var_cols
    long_df = (long_df.groupby(group_cols, observed=True, dropna=False)['TA_Lag1']
               .mean().reset_index())
    print(f"  Parcel TA averaged across runs: {len(long_df)} Subject x Session x Region rows")

    for var in var_cols:
        vals = sorted(long_df[var].dropna().astype(str).unique())
        long_df[var] = pd.Categorical(long_df[var], categories=vals)
    return long_df


def run_parcel_lmm(parcel_df, model_vars, ref_levels):
    if parcel_df.empty:
        return pd.DataFrame()
    parcel_vars   = [v for v in model_vars if v in parcel_df.columns]
    region_labels = sorted(parcel_df['Region_Label'].unique())
    print(f"  Running parcel-level LMM for {len(region_labels)} regions...")

    results = []
    for label in region_labels:
        sub = parcel_df[parcel_df['Region_Label'] == label].copy()
        sub = sub[['Subject'] + parcel_vars + ['TA_Lag1']].dropna()
        for var in [v for v in parcel_vars if v not in ('Session', 'Run')]:
            sub = sub[~sub[var].astype(str).isin(['NA', 'na', ''])]
        if sub.empty or sub['Subject'].nunique() < 2:
            continue

        formula, short_map, rename_map = build_formula(
            'TA_Lag1', parcel_vars, ref_levels)
        sub_model = sub.rename(columns={v: k for k, v in short_map.items()})
        try:
            result = smf.mixedlm(formula, data=sub_model,
                                 groups=sub_model['Subject']).fit(
                reml=True, method='powell', maxiter=500, disp=False)
            terms = rename_terms(list(result.params.index), rename_map)
            for term, est, se, z, p in zip(
                    terms, result.params.values, result.bse.values,
                    result.tvalues.values, result.pvalues.values):
                results.append({'Region_Label': label, 'term': term,
                                'estimate': est, 'std_error': se,
                                'z_value': z, 'p_value': p,
                                'n_obs': len(sub),
                                'n_subjects': sub['Subject'].nunique()})
        except Exception as e:
            results.append({'Region_Label': label,
                            'term': f'MODEL_FAILED: {e}',
                            'estimate': np.nan, 'std_error': np.nan,
                            'z_value': np.nan, 'p_value': np.nan,
                            'n_obs': len(sub),
                            'n_subjects': sub['Subject'].nunique()})

    return pd.DataFrame(results) if results else pd.DataFrame()


def plot_parcel_drivers(parcel_results, outdir, top_n=10):
    interaction_terms = [
        t for t in parcel_results['term'].dropna().unique()
        if ':' in str(t) and 'MODEL_FAILED' not in str(t)
    ]
    for term in interaction_terms:
        inter_df = (parcel_results[parcel_results['term'] == term]
                    .dropna(subset=['p_value'])
                    .sort_values('p_value').head(top_n))
        if inter_df.empty:
            continue
        safe_name = re.sub(r'[^\w]', '_', term)[:60]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh([str(int(r)) for r in inter_df['Region_Label']],
                -np.log10(inter_df['p_value'].clip(lower=1e-10)),
                color='steelblue', edgecolor='white')
        ax.axvline(-np.log10(0.05), color='red', linestyle='--',
                   linewidth=1.2, label='p=0.05')
        ax.set_xlabel(f'-log10(p-value)  [{term}]')
        ax.set_ylabel('Region Label')
        ax.set_title(f'Top {top_n} parcel drivers\n{term}', fontweight='bold')
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / f'parcel_drivers_{safe_name}.png',
                    dpi=300, bbox_inches='tight')
        fig.savefig(outdir / f'parcel_drivers_{safe_name}.pdf',
                    bbox_inches='tight')
        plt.close(fig)


# ============================================================
# README summary writer
# ============================================================

def write_summary(outdir, removed_df, mixed_results, model_vars,
                  ref_levels, region_label, outlier_scope):
    lines = [
        f"rs-fMRI analysis — {region_label}", "=" * 50, "",
        f"Model variables : {model_vars}",
        f"Reference levels: { {v: ref_levels.get(v) for v in model_vars} }",
        f"Formula         : metric ~ {'*'.join(model_vars)}  + (1|Subject)",
        f"Outlier scope   : {outlier_scope}",
        f"MAD multiplier  : {MAD_MULTIPLIER}",
        "", "Metrics analyzed:"
    ] + [f"  - {m}" for m in METRICS]

    lines += ["", "Outlier removal summary (run-metric values removed):"]
    if removed_df is None or removed_df.empty:
        lines.append("  - None")
    else:
        for metric in METRICS:
            sub_df = removed_df[removed_df['Metric'] == metric]
            if not sub_df.empty:
                lines.append(f"  {metric}: {len(sub_df)} removed")
                for _, r in sub_df.iterrows():
                    lines.append(
                        f"    {r['Subject']} | {r['Session']} | {r['Run']} | "
                        f"value={r['Value']:.4f}")

    sig = mixed_results[mixed_results['p_value'].notna() &
                        (mixed_results['p_value'] < 0.05)]
    lines += ["", "Terms with p < 0.05:"]
    lines += (["  - None"] if sig.empty else
              [f"  - {r['dataset']} | {r['metric']} | {r['term']} | "
               f"estimate={r['estimate']:.4f}, p={r['p_value']:.4g}"
               for _, r in sig.iterrows()])
    lines.append(f"\nOutputs written to: {outdir.resolve()}")
    (outdir / 'README_summary.txt').write_text("\n".join(lines))


# ============================================================
# Single-region-set analysis pipeline
# ============================================================

def run_analysis(metrics_csv, batch_dir, participants_tsv,
                 model_vars, outdir, label,
                 outlier_scope, filter_genotypes, n_dist_bins):
    print(f"\n{'='*60}")
    print(f"  ANALYSIS: {label}")
    print(f"  Model: metric ~ {'*'.join(model_vars)}  + (1|Subject)")
    print(f"  Outlier scope: {outlier_scope}")
    if filter_genotypes:
        print(f"  Genotype filter: {filter_genotypes}")
    print(f"{'='*60}")
    ensure_dir(outdir)

    # 1. Load + merge
    df = load_and_merge(metrics_csv, participants_tsv, model_vars)

    # 2. Genotype filter
    df = apply_genotype_filter(df, filter_genotypes)

    # Save full pre-outlier merged spreadsheet
    df.to_csv(outdir / 'all_subjects_merged_preoutlier.csv', index=False)
    print(f"  Saved pre-outlier merged data: "
          f"{df['Subject'].nunique()} subjects, {len(df)} rows")

    # 3. Prepare data (categoricals, reference levels)
    ref_levels = dict(REFERENCE_LEVELS)
    dat, ref_levels = prepare_data(df, model_vars, ref_levels, outdir)

    # 4. Outlier removal — per metric, independent
    removed_df = None
    if outlier_scope == 'within':
        dat, removed_df = remove_outliers_within(dat)
    elif outlier_scope == 'between':
        dat, removed_df = remove_outliers_between(dat)
    else:
        print("  No outlier removal applied (--outlier_scope not set)")

    if removed_df is not None and not removed_df.empty:
        removed_df.to_csv(outdir / 'outlier_runs_removed.csv', index=False)
        print(f"  Outlier log saved → outlier_runs_removed.csv")

    # Save run-level data after outlier removal (NaN where removed)
    dat.to_csv(outdir / 'all_subjects_runlevel_postoutlier.csv', index=False)

    # 5. Average runs within Subject x Session → LMM input
    dat_avg = average_runs(dat, model_vars)
    dat_avg.to_csv(outdir / 'all_subjects_averaged.csv', index=False)
    print(f"  Run-averaged data: "
          f"{dat_avg['Subject'].nunique()} subjects, {len(dat_avg)} rows")

    # 6. LMM on averaged data
    mixed_results = fit_all_models(dat_avg, model_vars, ref_levels,
                                   'run_averaged')
    mixed_results.to_csv(outdir / 'mixed_model_results.csv', index=False)

    interaction_summary = mixed_results[
        mixed_results['term'].astype(str).str.contains(':', na=False) &
        ~mixed_results['term'].astype(str).str.contains('MODEL_FAILED', na=False)
    ].copy()
    interaction_summary.to_csv(outdir / 'interaction_summary.csv', index=False)

    # 7. Plots — LMM trajectories
    for metric in METRICS:
        plot_metric(dat_avg, metric, model_vars,
                    f' — {label}',
                    outdir / f'{metric}_averaged')
        plot_raw_runlevel(dat, metric, model_vars,
                          outdir / f'{metric}_raw_runlevel')

    # 8. Distance-FC plot by genotype
    print(f"\n  Building distance-FC plot...")
    plot_distance_fc_by_group(batch_dir, dat_avg, outdir,
                               n_dist_bins, filter_genotypes)

    # 9. Parcel-level TA driver analysis
    print(f"\n  Running parcel-level TA driver analysis for {label}...")
    parcel_subdir = outdir / 'parcel_level_ta'
    ensure_dir(parcel_subdir)
    parcel_df = collect_parcel_ta(batch_dir, dat_avg, model_vars)
    if parcel_df.empty:
        print("  No parcel TA data found — skipping parcel analysis.")
    else:
        parcel_results = run_parcel_lmm(parcel_df, model_vars, ref_levels)
        if not parcel_results.empty:
            parcel_results.to_csv(
                parcel_subdir / 'parcel_lmm_all_terms.csv', index=False)
            parcel_results[
                parcel_results['term'].astype(str).str.contains(':', na=False) &
                ~parcel_results['term'].astype(str).str.contains(
                    'MODEL_FAILED', na=False)
            ].sort_values('p_value').to_csv(
                parcel_subdir / 'parcel_lmm_interaction_summary.csv',
                index=False)
            plot_parcel_drivers(parcel_results, parcel_subdir, top_n=10)

    # 10. README
    write_summary(outdir, removed_df, mixed_results, model_vars,
                  ref_levels, label, outlier_scope or 'none')
    print(f"\n  {label} analysis complete → {outdir.resolve()}")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='LMM analysis — SA/TA for one or more region sets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_lmm_isocortex_subcortex.py \\
      --analyses wholebrain \\
      --csv  /path/wholebrain.csv \\
      --batch /path/wb_batch \\
      --participants /path/participants.tsv \\
      --model_vars Genotype Session \\
      --outlier_scope within \\
      --filter_genotypes KI2APOE4 KI2APOE3 \\
      --n_dist_bins 20 \\
      --output_dir /path/output

  # Multiple analyses:
      --analyses isocortex subcortex wholebrain \\
      --csv  /path/iso.csv /path/sub.csv /path/wb.csv \\
      --batch /path/iso_batch /path/sub_batch /path/wb_batch
        """
    )
    parser.add_argument('--analyses', required=True, nargs='+',
                        help='Region sets to analyse (e.g. isocortex subcortex wholebrain). '
                             'Must match order of --csv and --batch.')
    parser.add_argument('--csv', required=True, nargs='+', dest='csvs',
                        help='Autocorrelation summary CSV per analysis (same order as --analyses).')
    parser.add_argument('--batch', required=True, nargs='+', dest='batches',
                        help='Batch output directory per analysis (same order as --analyses).')
    parser.add_argument('--participants', required=True,
                        help='TSV with Subject and any of: Treatment, Genotype, sex, session, run.')
    parser.add_argument('--model_vars', required=True, nargs='+',
                        help='Model variables — full factorial. '
                             'Recognised: Treatment, Genotype, Session. '
                             'Example: --model_vars Genotype Session')
    parser.add_argument('--outlier_scope', choices=['within', 'between'], default=None,
                        help="'within': per-subject MAD across runs (for subjects with many runs). "
                             "'between': cross-subject MAD on run-pair delta (for 2-run subjects). "
                             "Omit to skip outlier removal.")
    parser.add_argument('--filter_genotypes', nargs='+', default=None,
                        help='Only include these genotype(s). '
                             'Example: --filter_genotypes KI2APOE4 KI2APOE3')
    parser.add_argument('--n_dist_bins', type=int, default=20,
                        help='Number of distance bins for the distance-FC plot (default: 20).')
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()

    recognised = {'Treatment', 'Genotype', 'Session'}
    unknown    = [v for v in args.model_vars if v not in recognised]
    if unknown:
        parser.error(f"Unrecognised --model_vars: {unknown}. "
                     f"Recognised: {sorted(recognised)}")

    n = len(args.analyses)
    if len(args.csvs) != n:
        parser.error(f"--csv needs {n} path(s) to match --analyses, got {len(args.csvs)}")
    if len(args.batches) != n:
        parser.error(f"--batch needs {n} path(s) to match --analyses, got {len(args.batches)}")

    outdir = Path(args.output_dir)

    for analysis, csv_path, batch_path in zip(args.analyses, args.csvs, args.batches):
        run_analysis(
            metrics_csv      = csv_path,
            batch_dir        = batch_path,
            participants_tsv = args.participants,
            model_vars       = args.model_vars,
            outdir           = outdir / analysis,
            label            = analysis.capitalize(),
            outlier_scope    = args.outlier_scope,
            filter_genotypes = args.filter_genotypes or [],
            n_dist_bins      = args.n_dist_bins,
        )

    print(f"\nAll analyses complete. Outputs written to: {outdir.resolve()}\n")


if __name__ == '__main__':
    main()
