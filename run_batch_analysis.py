import os
import argparse
import pandas as pd
import glob
import re
import sys
from calculate_autocorrelation import load_and_validate_data, extract_centroids, calculate_ta, calculate_sa

def parse_bids_from_folder(folder_name):
    """
    Extracts Subject, Session, and Run from the specific folder naming convention:
    e.g., "_split_name_sub-AS160M1_dir-AP_ses-01_task-rest_run-2_part-mag_bold"
    """
    # Regex to find sub-..., ses-..., run-...
    # The convention seems to be key-value pairs separated by underscores or hyphens
    
    sub_match = re.search(r'(sub-[a-zA-Z0-9]+)', folder_name)
    ses_match = re.search(r'(ses-[a-zA-Z0-9]+)', folder_name)
    run_match = re.search(r'(run-[a-zA-Z0-9]+)', folder_name)
    
    sub = sub_match.group(1) if sub_match else "unknown"
    ses = ses_match.group(1) if ses_match else "unknown"
    run = run_match.group(1) if run_match else "unknown"
    
    return sub, ses, run

def find_matching_fc_file(base_fc_dir, sub, ses, run):
    """
    Searches for the corresponding FC matrix file.
    It ignores the 'dir-AP' vs 'dir-PA' difference by using wildcards or regex.
    """
    # Pattern to look for in the FC directory
    # Structure: .../analysis_main_wf/analysis_wf/_split_name_[ID_BLOCK]/FC_matrix/*.csv
    
    # We construct a regex pattern that matches the critical components
    # The folder name starts with _split_name_
    # Contains sub, ses, run
    # Might contain dir-AP or dir-PA (or other things)
    
    search_root = os.path.join(base_fc_dir, 'analysis_main_wf', 'analysis_wf')
    
    if not os.path.exists(search_root):
        print(f"  Warning: FC search root not found: {search_root}")
        return None
        
    # List all directories in the analysis root
    potential_dirs = [d for d in os.listdir(search_root) if os.path.isdir(os.path.join(search_root, d))]
    
    found_file = None
    
    for d in potential_dirs:
        # Check if this directory contains our subject, session, and run
        if sub in d and ses in d and run in d:
            # Found a candidate folder!
            # Now look for the CSV inside FC_matrix/
            fc_folder = os.path.join(search_root, d, 'FC_matrix')
            if os.path.exists(fc_folder):
                csv_files = glob.glob(os.path.join(fc_folder, '*.csv'))
                if csv_files:
                    found_file = csv_files[0] # Take the first one found
                    break
    
    return found_file

def main():
    parser = argparse.ArgumentParser(description="Batch Process TA and SA for Multiple Subjects")
    parser.add_argument("--root_dir", default=".", help="Root directory containing conf_correction_out and analysis_out_fcmatrix subdirectories (used as base path when --timeseries or --fc_matrix are not specified)")
    parser.add_argument("--atlas", required=True, help="Path to the Atlas NIfTI file")
    parser.add_argument("--output_dir", required=True, help="Directory to save aggregated results")
    parser.add_argument("--timeseries", default=None, help="Path to timeseries directory (overrides default path within root_dir)")
    parser.add_argument("--fc_matrix", default=None, help="Path to FC matrix directory (overrides default path within root_dir)")
    
    args = parser.parse_args()
    
    # Define Input Roots
    ts_root = args.timeseries or os.path.join(args.root_dir, 'conf_correction_out', 'confound_correction_datasink', 'cleaned_timeseries')
    fc_root = args.fc_matrix or os.path.join(args.root_dir, 'analysis_out_fcmatrix')
    
    if not os.path.exists(ts_root):
        print(f"Error: Timeseries directory not found at: {ts_root}")
        sys.exit(1)
    if not os.path.exists(fc_root):
        print(f"Error: FC Matrix directory not found at: {fc_root}")
        sys.exit(1)
        
    # Prepare Output
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    # List to store summary data
    summary_results = []
    
    # 1. Load Atlas (Once)
    print("Loading Atlas...")
    try:
        atlas_img = nib.load(args.atlas)
        atlas_data = atlas_img.get_fdata()
        atlas_affine = atlas_img.affine
    except Exception as e:
        print(f"Error loading atlas: {e}")
        sys.exit(1)
        
    # 2. Find Timeseries Files
    print(f"\nScanning for subjects in {ts_root}...")
    
    # Check if the directory contains CSV files directly (new behavior)
    direct_csv_files = glob.glob(os.path.join(ts_root, '*.csv'))
    
    if direct_csv_files:
        # New behavior: Process CSV files directly from the directory
        print(f"Found {len(direct_csv_files)} CSV files in directory.")
        ts_items = [(None, csv_file) for csv_file in direct_csv_files]
    else:
        # Original behavior: Look for nested folder structure
        # The structure is .../cleaned_timeseries/_split_name_.../*.csv
        ts_folders = [d for d in os.listdir(ts_root) if os.path.isdir(os.path.join(ts_root, d)) and d.startswith('_split_name_')]
        print(f"Found {len(ts_folders)} potential scans in nested folders.")
        ts_items = [(folder, None) for folder in ts_folders]
    
    for folder, direct_csv in ts_items:
        if direct_csv:
            # Extract BIDS info from the CSV filename
            csv_basename = os.path.basename(direct_csv)
            sub, ses, run = parse_bids_from_folder(csv_basename)
            print(f"\nProcessing: {sub} | {ses} | {run}")
            ts_path = direct_csv
        else:
            # Extract BIDS info from folder name
            sub, ses, run = parse_bids_from_folder(folder)
            print(f"\nProcessing: {sub} | {ses} | {run}")
            
            # Locate Timeseries CSV in folder
            ts_folder_path = os.path.join(ts_root, folder)
            ts_files = glob.glob(os.path.join(ts_folder_path, '*.csv'))
            
            if not ts_files:
                print("  Skipping: No CSV found in timeseries folder.")
                continue
                
            ts_path = ts_files[0]
        
        # Locate FC Matrix CSV
        fc_path = find_matching_fc_file(fc_root, sub, ses, run)
        
        if not fc_path:
            print(f"  Skipping: Could not find matching FC matrix for {sub} {ses} {run}")
            continue
            
        print(f"  TS: {os.path.basename(ts_path)}")
        print(f"  FC: {os.path.basename(fc_path)}")
        
        # Create output subdirectory for this scan
        scan_id = f"{sub}_{ses}_{run}"
        scan_out_dir = os.path.join(args.output_dir, scan_id)
        if not os.path.exists(scan_out_dir):
            os.makedirs(scan_out_dir)
            
        try:
            # --- RUN ANALYSIS ---
            # We call the functions directly instead of subprocess for speed/efficiency
            # 1. Validate
            # Note: We reuse the loaded atlas data to avoid reloading 100 times
            # However, `load_and_validate_data` currently loads the atlas inside it.
            # To avoid refactoring the main script too much, we will just call load_and_validate_data
            # but pass the path. It's slightly inefficient but safer to reuse tested code.
            # Actually, let's optimize slightly: if we look at `load_and_validate_data`
            # it returns atlas_data. 
            # We can just run it as is.
            
            fc_matrix, ts_data, region_labels, _, _ = load_and_validate_data(fc_path, ts_path, args.atlas)
            
            # 2. Extract Centroids (We can cache this technically, but it's fast)
            centroids = extract_centroids(region_labels, atlas_data, atlas_affine)
            
            # 3. Calculate TA
            global_ta, ta_df = calculate_ta(ts_data, region_labels)
            
            # Add coordinates
            ta_df['x'] = centroids[:, 0]
            ta_df['y'] = centroids[:, 1]
            ta_df['z'] = centroids[:, 2]
            
            # 4. Calculate SA
            sa_lambda, sa_inf, bin_x, bin_y, fit_y = calculate_sa(fc_matrix, centroids)
            
            # 5. Save Individual Results
            ta_df.to_csv(os.path.join(scan_out_dir, 'ta_results.csv'), index=False)
            
            if bin_x is not None:
                curve_df = pd.DataFrame({
                    'Distance_mm': bin_x,
                    'Mean_Correlation': bin_y,
                    'Fitted_Model': fit_y if fit_y is not None else [np.nan]*len(bin_x)
                })
                curve_df.to_csv(os.path.join(scan_out_dir, 'sa_curve_fit.csv'), index=False)
            
            # 6. Append to Summary
            summary_results.append({
                'Subject': sub,
                'Session': ses,
                'Run': run,
                'Global_TA': global_ta,
                'SA_Lambda': sa_lambda,
                'SA_Infinity': sa_inf
            })
            
            print("  ✅ Success")
            
        except Exception as e:
            print(f"  ❌ Error processing scan: {e}")
            import traceback
            traceback.print_exc()

    # --- FINAL SUMMARY ---
    if summary_results:
        print("\nWriting Master Summary...")
        summary_df = pd.DataFrame(summary_results)
        
        # Reorder columns nicely
        cols = ['Subject', 'Session', 'Run', 'Global_TA', 'SA_Lambda', 'SA_Infinity']
        summary_df = summary_df[cols]
        
        output_file = os.path.join(args.output_dir, 'all_subjects_autocorrelation.csv')
        summary_df.to_csv(output_file, index=False)
        print(f"Done! Master results saved to: {output_file}")
    else:
        print("\nNo scans were successfully processed.")

if __name__ == "__main__":
    # Needed for direct calls to functions in the other script
    import nibabel as nib
    import numpy as np
    main()
