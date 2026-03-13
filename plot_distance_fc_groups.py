"""
plot_distance_fc_groups.py

Plots binned functional connectivity vs. inter-region distance for multiple
groups overlaid on the same axes.

Input data comes from the sa_curve_fit.csv files produced by the batch
autocorrelation scripts (batch_autocorrelation_wholebrain.py etc.).
Each sa_curve_fit.csv contains pre-binned distance-FC data for one scan,
with columns: Distance_mm, Mean_Correlation, Fitted_Model.

Groups are defined by a participants TSV (the same one used in run_lmm_*).
Any column in that TSV can be used as the grouping variable (e.g. Genotype,
Treatment, sex).

Each group is plotted as:
  - Points at each bin centre showing the cross-subject mean FC
  - Vertical error bars showing ±1 SEM across subjects
  - An exponential decay fit line through the group mean points

Usage — single grouping variable:
    python plot_distance_fc_groups.py \
        --batch_dir  /path/to/batch_output/ \
        --participants /path/to/participants_mapping.tsv \
        --group_by   Genotype \
        --groups     KI2APOE3 KI2APOE4 \
        --output     distance_fc_genotype.png

Usage — filter to a subset (e.g. one session only):
    python plot_distance_fc_groups.py \
        --batch_dir  /path/to/batch_output/ \
        --participants /path/to/participants_mapping.tsv \
        --group_by   Treatment \
        --groups     Vehicle Drug \
        --filter_col Session \
        --filter_val ses-06m \
        --output     distance_fc_treatment_ses06m.png

Usage — cohort mode (batch_dir contains cohort subfolders):
    python plot_distance_fc_groups.py \
        --batch_dir  /path/to/batch_output/ \
        --participants /path/to/participants_mapping.tsv \
        --group_by   Genotype \
        --output     distance_fc_genotype.png \
        --cohort_mode
"""

import argparse
import os
import re
import glob
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Default colour palette — colour-blind friendly
# ---------------------------------------------------------------------------
DEFAULT_COLORS = [
    "#2271B2",   # blue
    "#E69F00",   # orange
    "#009E73",   # green
    "#CC79A7",   # pink
    "#56B4E9",   # sky blue
    "#D55E00",   # vermillion
    "#F0E442",   # yellow
    "#000000",   # black
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_bids_subject(path):
    """Extract sub- identifier from a file path."""
    m = re.search(r'(sub-[a-zA-Z0-9]+)', path)
    return m.group(1) if m else None


def load_participants(tsv_path):
    """Load participants TSV, strip whitespace from column names and values."""
    df = pd.read_csv(tsv_path, sep='\t')
    df.columns = df.columns.str.strip()
    for col in df.select_dtypes(include='object').columns:
        df[col] = df[col].str.strip()
    # Normalise subject column name
    for candidate in ['Subject', 'subject', 'participant_id', 'Participant']:
        if candidate in df.columns:
            df = df.rename(columns={candidate: 'Subject'})
            break
    return df


def collect_sa_curves(batch_dir, cohort_mode=False):
    """
    Recursively finds all sa_curve_fit.csv files under batch_dir.
    Returns a list of dicts: {subject, path}
    """
    records = []
    pattern = os.path.join(batch_dir, '**', 'sa_curve_fit.csv')
    for path in glob.glob(pattern, recursive=True):
        sub = parse_bids_subject(path)
        if sub:
            records.append({'Subject': sub, 'path': path})
    return records


def load_and_bin_curves(records, n_bins=20):
    """
    Loads all sa_curve_fit.csv files and resamples them onto a common
    distance grid using n_bins equal-width bins across the global distance
    range. Returns a DataFrame with columns:
        Subject, bin_center, mean_fc
    """
    # First pass — find global distance range
    all_dists = []
    loaded = []
    for rec in records:
        try:
            df = pd.read_csv(rec['path'])
            df = df.dropna(subset=['Distance_mm', 'Mean_Correlation'])
            if len(df) < 2:
                continue
            all_dists.extend(df['Distance_mm'].tolist())
            loaded.append((rec['Subject'], df))
        except Exception as e:
            print(f"  Warning: could not load {rec['path']}: {e}")

    if not loaded:
        raise ValueError("No valid sa_curve_fit.csv files could be loaded.")

    d_min = np.min(all_dists)
    d_max = np.max(all_dists)
    bin_edges = np.linspace(d_min, d_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    rows = []
    for sub, df in loaded:
        dists = df['Distance_mm'].values
        fcs   = df['Mean_Correlation'].values
        for i, center in enumerate(bin_centers):
            mask = (dists >= bin_edges[i]) & (dists < bin_edges[i + 1])
            vals = fcs[mask]
            vals = vals[~np.isnan(vals)]
            if len(vals) > 0:
                rows.append({
                    'Subject':    sub,
                    'bin_center': center,
                    'mean_fc':    float(np.mean(vals))
                })

    return pd.DataFrame(rows), bin_centers


def exp_decay(x, a, lam, c):
    """Exponential decay: a * exp(-x / lam) + c"""
    return a * np.exp(-x / lam) + c


def fit_exp_decay(bin_centers, group_means):
    """Fit exponential decay to group mean curve. Returns fitted y values."""
    try:
        a0   = group_means[0] - group_means[-1]
        lam0 = (bin_centers[-1] - bin_centers[0]) / 3.0
        c0   = group_means[-1]
        popt, _ = curve_fit(
            exp_decay, bin_centers, group_means,
            p0=[a0, lam0, c0],
            bounds=([-2, 0.001, -2], [2, 500, 2]),
            maxfev=10000
        )
        x_fine = np.linspace(bin_centers[0], bin_centers[-1], 300)
        return x_fine, exp_decay(x_fine, *popt)
    except Exception:
        # Fall back to linear
        m, b = np.polyfit(bin_centers, group_means, 1)
        x_fine = np.linspace(bin_centers[0], bin_centers[-1], 300)
        return x_fine, m * x_fine + b


# ---------------------------------------------------------------------------
# Main plot function
# ---------------------------------------------------------------------------

def plot_groups(binned_df, bin_centers, group_assignments, group_order,
                group_by, colors, output_path, title, xlabel, ylabel,
                figsize, point_size, show_fit, show_sem, filter_desc):
    """
    binned_df      : DataFrame with Subject, bin_center, mean_fc
    bin_centers    : array of bin centre distances
    group_assignments : dict mapping Subject -> group label
    group_order    : ordered list of group labels to plot
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Offset bin centers slightly per group so points don't overlap
    n_groups = len(group_order)
    offsets  = np.linspace(-0.3, 0.3, n_groups) if n_groups > 1 else [0]

    for gi, group_label in enumerate(group_order):
        color   = colors[gi % len(colors)]
        offset  = offsets[gi]

        # Subjects in this group
        subjects_in_group = [s for s, g in group_assignments.items()
                             if g == group_label]
        if not subjects_in_group:
            print(f"  Warning: no subjects found for group '{group_label}'")
            continue

        group_data = binned_df[binned_df['Subject'].isin(subjects_in_group)]
        if group_data.empty:
            print(f"  Warning: no curve data found for group '{group_label}'")
            continue

        print(f"  Group '{group_label}': {len(subjects_in_group)} subjects, "
              f"{group_data['Subject'].nunique()} with curve data")

        # Aggregate across subjects per bin
        agg = (group_data
               .groupby('bin_center', sort=True)['mean_fc']
               .agg(['mean', 'std', 'count'])
               .reset_index())
        agg['sem'] = agg['std'] / np.sqrt(agg['count'])
        agg = agg.dropna(subset=['mean'])

        x_plot   = agg['bin_center'].values + offset
        y_plot   = agg['mean'].values
        sem_plot = agg['sem'].values

        # Points
        ax.scatter(
            x_plot, y_plot,
            s=point_size,
            color=color,
            zorder=3,
            label=f"{group_label} (n={len(subjects_in_group)})"
        )

        # Error bars
        if show_sem:
            ax.errorbar(
                x_plot, y_plot,
                yerr=sem_plot,
                fmt='none',
                color=color,
                capsize=3,
                linewidth=1.2,
                zorder=2
            )

        # Exponential fit line
        if show_fit and len(agg) >= 3:
            x_fine, y_fine = fit_exp_decay(agg['bin_center'].values, y_plot)
            ax.plot(x_fine, y_fine,
                    color=color, linewidth=2.0, zorder=2,
                    alpha=0.85)

    ax.axhline(0, color='#999999', linewidth=0.8, linestyle='--', zorder=1)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    full_title = title
    if filter_desc:
        full_title += f'\n{filter_desc}'
    ax.set_title(full_title, fontsize=13, fontweight='bold')

    ax.legend(fontsize=10, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()

    # Save both PNG and PDF
    base = os.path.splitext(output_path)[0]
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    print(f"Plot saved to: {base}.pdf")

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot binned distance-FC curves for multiple groups overlaid"
    )
    parser.add_argument('--batch_dir', required=True,
                        help='Batch output directory containing sa_curve_fit.csv files')
    parser.add_argument('--participants', required=True,
                        help='Path to participants TSV file')
    parser.add_argument('--group_by', required=True,
                        help='Column in participants TSV to use for grouping '
                             '(e.g. Genotype, Treatment, sex)')
    parser.add_argument('--groups', nargs='+', default=None,
                        help='Ordered list of group values to include and plot. '
                             'If not supplied, all unique values are used.')
    parser.add_argument('--filter_col', default=None,
                        help='Optional: column to filter participants on '
                             '(e.g. Session)')
    parser.add_argument('--filter_val', default=None,
                        help='Value to keep for --filter_col (e.g. ses-06m)')
    parser.add_argument('--output', default='distance_fc_groups.png',
                        help='Output plot path (PNG; PDF also saved automatically)')
    parser.add_argument('--n_bins', type=int, default=20,
                        help='Number of distance bins (default: 20)')
    parser.add_argument('--title', default='Distance-Dependent Functional Connectivity',
                        help='Plot title')
    parser.add_argument('--xlabel', default='Inter-region Distance (mm)',
                        help='X-axis label')
    parser.add_argument('--ylabel', default='Mean Functional Connectivity (r)',
                        help='Y-axis label')
    parser.add_argument('--figsize', nargs=2, type=float, default=[8, 5],
                        metavar=('WIDTH', 'HEIGHT'),
                        help='Figure size in inches (default: 8 5)')
    parser.add_argument('--point_size', type=float, default=40,
                        help='Scatter point size (default: 40)')
    parser.add_argument('--no_fit', action='store_true',
                        help='Disable exponential fit line overlay')
    parser.add_argument('--no_sem', action='store_true',
                        help='Disable SEM error bars')
    parser.add_argument('--colors', nargs='+', default=None,
                        help='Custom hex or named colours per group, '
                             'in the same order as --groups')
    parser.add_argument('--cohort_mode', action='store_true',
                        help='Batch dir contains cohort subfolders '
                             '(cohort_100 etc.) — recurse into them')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    # --- Load participants ---
    print("Loading participants TSV...")
    participants = load_participants(args.participants)
    print(f"  {len(participants)} rows, columns: {participants.columns.tolist()}")

    if 'Subject' not in participants.columns:
        raise ValueError("Could not find a Subject column in the participants TSV.")

    # Case-insensitive column matching for group_by
    col_map = {c.lower(): c for c in participants.columns}
    if args.group_by.lower() not in col_map:
        raise ValueError(f"Column '{args.group_by}' not found in participants TSV. "
                         f"Available: {participants.columns.tolist()}")
    args.group_by = col_map[args.group_by.lower()]

    # --- Optional row filter ---
    filter_desc = None
    if args.filter_col and args.filter_val:
        if args.filter_col.lower() not in col_map:
            raise ValueError(f"--filter_col '{args.filter_col}' not in participants TSV.")
        args.filter_col = col_map[args.filter_col.lower()]
        participants = participants[
            participants[args.filter_col].astype(str) == args.filter_val
        ].copy()
        filter_desc = f"{args.filter_col} = {args.filter_val}"
        print(f"  After filter ({filter_desc}): {len(participants)} rows")

    # --- Build group assignments (case-insensitive value matching) ---
    # Build a map from lowercase TSV value -> original TSV value
    tsv_values = participants[args.group_by].dropna().unique().tolist()
    tsv_val_map = {str(v).lower(): str(v) for v in tsv_values}

    if args.groups:
        # Map user-supplied group names to actual TSV values
        resolved_groups = []
        for g in args.groups:
            if str(g).lower() in tsv_val_map:
                resolved_groups.append(tsv_val_map[str(g).lower()])
            else:
                print(f"  Warning: group '{g}' not found in TSV column "
                      f"'{args.group_by}'. Available: {list(tsv_val_map.values())}")
        group_order = resolved_groups
    else:
        group_order = sorted(tsv_values)

    participants = participants[
        participants[args.group_by].isin(group_order)
    ].copy()

    group_assignments = dict(
        zip(participants['Subject'].astype(str),
            participants[args.group_by].astype(str))
    )
    print(f"  Groups: {group_order}")
    for g in group_order:
        n = sum(1 for v in group_assignments.values() if v == g)
        print(f"    {g}: {n} subjects")

    # --- Collect sa_curve_fit.csv files ---
    print(f"\nSearching for sa_curve_fit.csv files in: {args.batch_dir}")
    records = collect_sa_curves(args.batch_dir, cohort_mode=args.cohort_mode)
    print(f"  Found {len(records)} sa_curve_fit.csv files")

    if not records:
        raise ValueError(
            f"No sa_curve_fit.csv files found under {args.batch_dir}. "
            "Check that the batch script has been run and --batch_dir is correct."
        )

    # Report which subjects have curve data
    found_subs = set(r['Subject'] for r in records)
    assigned_subs = set(group_assignments.keys())
    matched = found_subs & assigned_subs
    print(f"  Subjects with curve data AND group assignment: {len(matched)}")
    missing_data = assigned_subs - found_subs
    if missing_data:
        print(f"  Subjects in participants TSV but no curve data: "
              f"{sorted(missing_data)}")

    # --- Load and bin curves ---
    print(f"\nLoading and binning curves ({args.n_bins} bins)...")
    # Filter records to only subjects with a group assignment
    records = [r for r in records if r['Subject'] in group_assignments]
    binned_df, bin_centers = load_and_bin_curves(records, n_bins=args.n_bins)
    print(f"  Binned data: {len(binned_df)} rows across "
          f"{binned_df['Subject'].nunique()} subjects")

    # --- Colours ---
    if args.colors:
        colors = args.colors
    else:
        colors = DEFAULT_COLORS

    # --- Plot ---
    print(f"\nPlotting {len(group_order)} groups...")
    plot_groups(
        binned_df=binned_df,
        bin_centers=bin_centers,
        group_assignments=group_assignments,
        group_order=group_order,
        group_by=args.group_by,
        colors=colors,
        output_path=args.output,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        figsize=tuple(args.figsize),
        point_size=args.point_size,
        show_fit=not args.no_fit,
        show_sem=not args.no_sem,
        filter_desc=filter_desc,
    )


if __name__ == '__main__':
    main()