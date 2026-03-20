# LMM SA/TA Analysis Pipeline

Linear mixed model analysis of spatial autocorrelation (SA) and temporal
autocorrelation (TA) metrics from resting-state fMRI data. Supports isocortex,
subcortex, and whole-brain region sets, with optional genotype filtering,
per-metric outlier removal, run averaging, and distance-dependent FC plots.

---

## Requirements

```
numpy
pandas
matplotlib
scipy
statsmodels
```

---

## Quick Start

```bash
python run_lmm_isocortex_subcortex.py \
    --analyses isocortex subcortex \
    --csv  /path/to/isocortex/all_subjects_isocortex_autocorrelation.csv \
           /path/to/subcortex/all_subjects_subcortex_autocorrelation.csv \
    --batch /path/to/isocortex/ \
            /path/to/subcortex/ \
    --participants /path/to/participants_treatment_map.tsv \
    --model_vars Treatment Session \
    --outlier_scope between \
    --n_dist_bins 20 \
    --output_dir /path/to/lmm_outputs/
```

---

## All Command-Line Flags

| Flag | Required | Description |
|------|----------|-------------|
| `--analyses` | Yes | One or more region set names (e.g. `isocortex subcortex wholebrain`). Each must have a matching `--csv` and `--batch` in the same order. |
| `--csv` | Yes | Path to the autocorrelation summary CSV for each analysis, in the same order as `--analyses`. |
| `--batch` | Yes | Path to the batch output directory for each analysis, in the same order as `--analyses`. Used to find `ta_results.csv` and `sa_curve_fit.csv` files for parcel-level and distance-FC analyses. |
| `--participants` | Yes | Tab-separated file with at minimum a `Subject` column, plus any of: `Treatment`, `Genotype`, `sex`, `session`, `run`. Column names are case-insensitive. |
| `--model_vars` | Yes | Variables for the fixed effects. Full factorial (all main effects + all interactions) is always used. Recognised values: `Treatment`, `Genotype`, `Session`. Example: `--model_vars Treatment Session` |
| `--outlier_scope` | No | Either `within` or `between`. Omit to skip outlier removal entirely. See details below. |
| `--filter_genotypes` | No | Restrict the analysis to specific genotype(s). Example: `--filter_genotypes KI2APOE4 KI2APOE3` |
| `--n_dist_bins` | No | Number of distance bins for the distance-dependent FC plot. Default: `20`. |
| `--output_dir` | Yes | Root output directory. Each analysis gets its own subdirectory inside it. |

---

## Outlier Removal

Outliers are identified and removed **independently per metric** — a run flagged
as an outlier for `Global_TA` is only removed from `Global_TA`; its `SA_Lambda`
and `SA_Infinity` values are kept. The MAD multiplier is 3.0 (no normalisation
factor applied).

### `--outlier_scope within`
Designed for subjects with **many runs** (e.g. 4 runs per session).

For each Subject × Session:
1. Compute the median and MAD of that subject's run values for that metric.
2. Flag any individual run whose value falls outside `median ± 3 × MAD`.
3. Set that run's value to NaN for that metric only.

### `--outlier_scope between`
Designed for subjects with **2 runs per session**.

For each Subject × Session:
1. Compute the absolute difference (range) between the 2 runs per metric.
2. Collect all ranges across all subjects.
3. Flag any Subject × Session where range > `median(ranges) + 3 × MAD(ranges)`.
4. For flagged sessions, remove the single run furthest from that subject's mean.

All removed values are logged to `outlier_runs_removed.csv` in each analysis
output folder.

---

## Run Averaging

After outlier removal (or immediately after loading if no outlier removal is
requested), all metric values are **averaged across runs within each
Subject × Session**. This averaged data is what is passed into the linear mixed
model. The run-level data (with outlier NaNs) is saved separately as
`all_subjects_runlevel_postoutlier.csv`.

---

## Linear Mixed Model

The model is fit using `statsmodels` REML with a Powell optimiser:

```
metric ~ V1 * V2 * ... + (1 | Subject)
```

Where `V1`, `V2`, etc. are the variables specified by `--model_vars`. Full
factorial means all main effects and all pairwise (and higher-order)
interactions are included. One model is fit per metric (`Global_TA`,
`SA_Lambda`, `SA_Infinity`).

Reference levels:
- `Treatment`: Vehicle
- `Genotype`: first alphabetically (if not overridden)
- `Session`: earliest session (sorted numerically)

---

## Distance-Dependent FC Plot

After the LMM, the script reads `sa_curve_fit.csv` files from the batch
directory, bins inter-region distances into `--n_dist_bins` equal-width bins,
and plots the mean FC ± SEM per bin as points, with an exponential decay fit
line overlaid. All genotype groups appear on the same plot.

Output: `distance_fc_by_genotype.png` / `.pdf`

If no `sa_curve_fit.csv` files are found in the batch directory, this step is
skipped with a warning.

---

## Parcel-Level TA Driver Analysis

After the global LMM, a separate per-region LMM is run for `TA_Lag1` using
`ta_results.csv` files found recursively in the batch directory. For each
region, the same model formula is fit, and interaction terms are sorted by
p-value. The top 10 regions per interaction term are plotted as bar charts
(`-log10(p)` on the x-axis).

Outputs saved to `parcel_level_ta/` subdirectory:
- `parcel_lmm_all_terms.csv`
- `parcel_lmm_interaction_summary.csv`
- `parcel_drivers_<term>.png/.pdf`

---

## Output Files

For each analysis (e.g. `lmm_outputs/isocortex/`):

| File | Description |
|------|-------------|
| `all_subjects_merged_preoutlier.csv` | Fully merged data before any outlier removal |
| `all_subjects_runlevel_postoutlier.csv` | Run-level data after outlier NaNs applied |
| `outlier_runs_removed.csv` | Log of every removed run-metric value with reason |
| `all_subjects_averaged.csv` | Run-averaged data — primary LMM input |
| `subject_summary.txt` | Subject count and group breakdown |
| `mixed_model_results.csv` | Full LMM results for all metrics and terms |
| `interaction_summary.csv` | Interaction terms only, sorted by p-value |
| `{metric}_averaged.png/.pdf` | Trajectory plot (mean ± SEM per group per session) |
| `{metric}_raw_runlevel.png/.pdf` | Run-level scatter with removed points flagged |
| `distance_fc_by_genotype.png/.pdf` | Distance-binned FC plot by genotype group |
| `README_summary.txt` | Auto-generated run summary with flagged outliers and significant terms |
| `parcel_level_ta/` | Parcel-level TA LMM results and driver plots |

---

## Participants TSV Format

The TSV must have a `Subject` column. Additional columns are optional but
recognised:

```
Subject	Treatment	Genotype	sex	session	run
sub-01	Drug	KI2APOE4	M	ses-06m	run-1
sub-01	Drug	KI2APOE4	M	ses-06m	run-2
sub-02	Vehicle	KI2APOE3	F	ses-06m	run-1
```

- Column names are matched case-insensitively.
- If `session` and `run` columns are present, merging is done on
  Subject + Session + Run to avoid row duplication.
- If `Treatment` or `Genotype` are absent, those columns are filled with
  blank and excluded from the model.

---

## Example Commands

### Lecanemab project — Treatment × Session, 2 runs, between-subject outlier removal
```bash
python run_lmm_isocortex_subcortex.py \
    --analyses isocortex subcortex \
    --csv  /path/isocortex/all_subjects_isocortex_autocorrelation.csv \
           /path/subcortex/all_subjects_subcortex_autocorrelation.csv \
    --batch /path/isocortex/ \
            /path/subcortex/ \
    --participants /path/participants_treatment_map.tsv \
    --model_vars Treatment Session \
    --outlier_scope between \
    --output_dir /path/lmm_outputs/
```

### Aya project — Genotype × Session, 4 runs, within-subject outlier removal, genotype filter
```bash
python run_lmm_isocortex_subcortex.py \
    --analyses wholebrain \
    --csv  /path/wholebrain_analysis/all_subjects_wholebrain_autocorrelation.csv \
    --batch /path/wholebrain_analysis/ \
    --participants /path/participants_mapping.tsv \
    --model_vars Genotype Session \
    --outlier_scope within \
    --filter_genotypes KI2APOE4 KI2APOE3 \
    --n_dist_bins 20 \
    --output_dir /path/wholebrain_analysis/lmm_model/
```

### Three region sets at once
```bash
python run_lmm_isocortex_subcortex.py \
    --analyses isocortex subcortex wholebrain \
    --csv  /path/iso.csv /path/sub.csv /path/wb.csv \
    --batch /path/iso/ /path/sub/ /path/wb/ \
    --participants /path/participants.tsv \
    --model_vars Treatment Session \
    --outlier_scope between \
    --output_dir /path/lmm_outputs/
```
