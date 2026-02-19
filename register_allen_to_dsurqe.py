import argparse
import subprocess
import sys
import os


def run_registration(fixed_image, moving_image, output_prefix):
    """
    Runs ANTs registration (Rigid + Affine + SyN) to register the Allen Brain Atlas
    (moving) to the DSURQE atlas template (fixed).

    Parameters
    ----------
    fixed_image : str
        Path to the fixed image (DSURQE_40micron_average.nii.gz).
    moving_image : str
        Path to the moving image (zoaverage_template_50.nii).
    output_prefix : str
        Output prefix for registration files (e.g. 'full_registration_').
    """
    warped_output = output_prefix + "Warped.nii.gz"

    cmd = [
        "antsRegistration",
        "--dimensionality", "3",
        "--float", "0",
        "--output", f"[{output_prefix},{warped_output}]",
        "--interpolation", "Linear",
        "--use-histogram-matching", "0",
        "--winsorize-image-intensities", "[0.005,0.995]",
        "--initial-moving-transform", f"[{fixed_image},{moving_image},1]",
        # Rigid stage
        "--transform", "Rigid[0.1]",
        "--metric", f"MI[{fixed_image},{moving_image},0.25]",
        "--convergence", "[1000x500x250x100,1e-6,10]",
        "--shrink-factors", "8x4x2x1",
        "--smoothing-sigmas", "3x2x1x0vox",
        # Affine stage
        "--transform", "Affine[0.1]",
        "--metric", f"CC[{fixed_image},{moving_image},1,4]",
        "--convergence", "[1000x500x250x100,1e-6,10]",
        "--shrink-factors", "8x4x2x1",
        "--smoothing-sigmas", "3x2x1x0vox",
        # SyN stage
        "--transform", "SyN[0.1,3,0]",
        "--metric", f"CC[{fixed_image},{moving_image},1,4]",
        "--convergence", "[100x70x50x20,1e-6,10]",
        "--shrink-factors", "8x4x2x1",
        "--smoothing-sigmas", "3x2x1x0vox",
        "--verbose", "1",
    ]

    print("Running antsRegistration command:")
    print(" ".join(cmd))
    print()

    result = subprocess.run(cmd, check=False)

    if result.returncode != 0:
        print(
            f"\nERROR: antsRegistration exited with code {result.returncode}. "
            "Check the output above for details."
        )
        sys.exit(result.returncode)

    print(f"\nRegistration complete. Warped image saved to: {warped_output}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Register the Allen Brain Atlas to the DSURQE atlas template "
            "using ANTs (Rigid + Affine + SyN)."
        )
    )
    parser.add_argument(
        "--fixed",
        default="DSURQE_40micron_average.nii.gz",
        help="Path to the fixed image (default: DSURQE_40micron_average.nii.gz).",
    )
    parser.add_argument(
        "--moving",
        default="zoaverage_template_50.nii",
        help="Path to the moving image (default: zoaverage_template_50.nii).",
    )
    parser.add_argument(
        "--output_prefix",
        default="full_registration_",
        help="Output prefix for registration files (default: full_registration_).",
    )

    args = parser.parse_args()

    for path, label in [(args.fixed, "Fixed image"), (args.moving, "Moving image")]:
        if not os.path.isfile(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    run_registration(args.fixed, args.moving, args.output_prefix)


if __name__ == "__main__":
    main()
