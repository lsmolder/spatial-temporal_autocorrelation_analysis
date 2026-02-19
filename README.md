# Autocorrelation Analysis Script

This script calculates Temporal Autocorrelation (TA) and Spatial Autocorrelation (SA) for brain imaging data, following the methods described in Shinn et al. (2023).

## Requirements

- Python 3
- See `requirements.txt` for dependencies.

## Usage

### Single Subject Analysis

To analyze a single subject's data:

```bash
python calculate_autocorrelation.py --fc_matrix <path_to_fc_matrix.csv> --timeseries <path_to_timeseries.csv> --atlas <path_to_atlas.nii.gz> --output <output_directory>
```

#### Arguments
- `--fc_matrix`: Path to the Functional Connectivity matrix CSV.
    - Format: First row and first column should be region labels.
- `--timeseries`: Path to the fMRI time series CSV.
    - Format: First row should be region labels. Data starts from the second row (no time column).
- `--atlas`: Path to the NIfTI atlas file (`.nii` or `.nii.gz`) defining the brain regions.
- `--output`: Directory where results will be saved.

### Batch Analysis

To analyze multiple subjects automatically, use the batch processing script. This script expects a specific directory structure common in BIDS-like preprocessing pipelines.

```bash
python run_batch_analysis.py --root_dir <path_to_data_root> --atlas <path_to_atlas.nii.gz> --output_dir <path_to_results>
```

#### Directory Structure Assumption
The script looks for data in the following locations relative to `--root_dir`:
1.  **Timeseries:** `<root_dir>/conf_correction_out/confound_correction_datasink/cleaned_timeseries/`
2.  **FC Matrices:** `<root_dir>/analysis_out_fcmatrix/analysis_main_wf/analysis_wf/`

It automatically pairs files based on Subject, Session, and Run IDs (e.g., `sub-01`, `ses-01`, `run-1`), handling variations in folder naming (like `dir-AP` vs `dir-PA`).

**Flexible Timeseries Input:**
The `--timeseries` flag now supports two directory structures:
- **Direct CSV files:** Point to a directory containing CSV files directly (e.g., `sub-01_ses-01_run-1.csv`). The script will extract BIDS information from the filenames.
- **Nested folders (legacy):** Point to a directory containing `_split_name_*` subdirectories, each with a CSV file inside.

#### Arguments
- `--root_dir`: The root folder containing the `conf_correction_out` and `analysis_out_fcmatrix` directories. This is used as the base path for finding data when `--timeseries` or `--fc_matrix` are not explicitly specified. Defaults to the current directory (`.`).
- `--atlas`: Path to the **common** NIfTI atlas file used for all subjects.
- `--output_dir`: Directory where the master summary and individual results will be saved.
- `--timeseries`: (Optional) Path to timeseries directory. If not provided, defaults to `<root_dir>/conf_correction_out/confound_correction_datasink/cleaned_timeseries/`.
- `--fc_matrix`: (Optional) Path to FC matrix directory. If not provided, defaults to `<root_dir>/analysis_out_fcmatrix`.

#### Batch Outputs
1.  **`all_subjects_autocorrelation.csv`**: A master CSV containing results for all scans.
    - Columns: `Subject`, `Session`, `Run`, `Global_TA`, `SA_Lambda`, `SA_Infinity`.
2.  **Subdirectories**: Inside the output directory, a folder is created for each scan (e.g., `sub-01_ses-01_run-1/`) containing the detailed per-region TA and curve fit data.

---

## Atlas Registration: Allen Brain Atlas → DSURQE

`register_allen_to_dsurqe.py` registers the scalable Allen Brain Atlas (`zoaverage_template_50.nii`) to the DSURQE atlas template (`DSURQE_40micron_average.nii.gz`) used in rabies preprocessing, using ANTs with a three-stage strategy: **Rigid → Affine → SyN**.

### Requirements

- [ANTs](https://github.com/ANTsX/ANTs) must be installed and `antsRegistration` available on your `PATH` (or `$PATH` on Unix systems).

### Usage

```bash
python register_allen_to_dsurqe.py \
    --fixed DSURQE_40micron_average.nii.gz \
    --moving zoaverage_template_50.nii \
    --output_prefix full_registration_
```

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `--fixed` | `DSURQE_40micron_average.nii.gz` | Path to the fixed image (DSURQE template). |
| `--moving` | `zoaverage_template_50.nii` | Path to the moving image (Allen Brain Atlas template). |
| `--output_prefix` | `full_registration_` | Prefix for all output transform and warped image files. |

### Outputs

- `<output_prefix>Warped.nii.gz` – Allen Brain Atlas warped into DSURQE space.
- `<output_prefix>0GenericAffine.mat` – Rigid + Affine transform.
- `<output_prefix>1Warp.nii.gz` – SyN forward warp field.
- `<output_prefix>1InverseWarp.nii.gz` – SyN inverse warp field.

### Registration Strategy

| Stage | Transform | Metric | Convergence | Shrink Factors | Smoothing |
|---|---|---|---|---|---|
| 1 | Rigid[0.1] | MI (weight 0.25) | 1000×500×250×100, tol 1e-6 | 8×4×2×1 | 3×2×1×0 vox |
| 2 | Affine[0.1] | CC (radius 4) | 1000×500×250×100, tol 1e-6 | 8×4×2×1 | 3×2×1×0 vox |
| 3 | SyN[0.1,3,0] | CC (radius 4) | 100×70×50×20, tol 1e-6 | 8×4×2×1 | 3×2×1×0 vox |

---

## Methodology

### Temporal Autocorrelation (TA)
Calculated as the lag-1 autocorrelation of the time series for each region. Global TA is the mean of all regional TAs.

### Spatial Autocorrelation (SA)
Calculated by fitting an exponential decay model to the relationship between Functional Connectivity (correlation) and Euclidean distance between regions.

Model: $y = SA_{inf} + (1 - SA_{inf}) \times e^{-x / SA_{\lambda}}$

- $SA_{\lambda}$: Length scale (rate of decay).
- $SA_{inf}$: Asymptote (baseline correlation at long distances).
