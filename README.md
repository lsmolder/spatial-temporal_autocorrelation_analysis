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

#### Arguments
- `--root_dir`: The root folder containing the `conf_correction_out` and `analysis_out_fcmatrix` directories. Defaults to the current directory (`.`).
- `--atlas`: Path to the **common** NIfTI atlas file used for all subjects.
- `--output_dir`: Directory where the master summary and individual results will be saved.

#### Batch Outputs
1.  **`all_subjects_autocorrelation.csv`**: A master CSV containing results for all scans.
    - Columns: `Subject`, `Session`, `Run`, `Global_TA`, `SA_Lambda`, `SA_Infinity`.
2.  **Subdirectories**: Inside the output directory, a folder is created for each scan (e.g., `sub-01_ses-01_run-1/`) containing the detailed per-region TA and curve fit data.

## Methodology

### Temporal Autocorrelation (TA)
Calculated as the lag-1 autocorrelation of the time series for each region. Global TA is the mean of all regional TAs.

### Spatial Autocorrelation (SA)
Calculated by fitting an exponential decay model to the relationship between Functional Connectivity (correlation) and Euclidean distance between regions.

Model: $y = SA_{inf} + (1 - SA_{inf}) \times e^{-x / SA_{\lambda}}$

- $SA_{\lambda}$: Length scale (rate of decay).
- $SA_{inf}$: Asymptote (baseline correlation at long distances).
