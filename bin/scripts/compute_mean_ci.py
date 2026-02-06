import pandas as pd
import numpy as np
import argparse
import os
import sys

ALPHA = 0.05
RANDOM_SEED = 0

def median_ci(series, n_boot, alpha=ALPHA, seed=RANDOM_SEED):
    data = series.dropna().values
    if len(data) == 0:
        return np.nan, np.nan, np.nan

    med = np.median(data)
    rng = np.random.default_rng(seed)

    boot_meds = []
    for _ in range(n_boot):
        sample = rng.choice(data, size=len(data), replace=True)
        boot_meds.append(np.median(sample))

    lower = np.percentile(boot_meds, 100 * alpha / 2.0)
    upper = np.percentile(boot_meds, 100 * (1.0 - alpha / 2.0))

    return med, lower, upper

def main():
    parser = argparse.ArgumentParser(
        description="Compute median and bootstrap CI for each numeric column."
    )
    parser.add_argument("input_csv", help="Path to input CSV file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    parser.add_argument(
        "--n-boot",
        type=int,
        default=2000,
        help="Number of bootstrap resamples (positive integer, default: 2000)"
    )
    args = parser.parse_args()

    # Guard: check input file exists and is a file
    if not os.path.isfile(args.input_csv):
        print(f"Error: input file '{args.input_csv}' does not exist or is not a file.", file=sys.stderr)
        sys.exit(1)

    # Guard: check n_boot is positive
    if args.n_boot <= 0:
        print("Error: --n-boot must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    input_csv = args.input_csv
    output_csv = args.output_csv
    n_boot = args.n_boot

    # Load data
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"Error: failed to read '{input_csv}': {e}", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("Error: input CSV has no rows.", file=sys.stderr)
        sys.exit(1)

    # Treat 'sample' as index if present
    if "sample" in df.columns:
        df = df.set_index("sample")

    results = []

    # Loop over numeric columns
    numeric_cols = [col for col in df.columns if np.issubdtype(df[col].dtype, np.number)]
    if not numeric_cols:
        print("Error: input CSV has no numeric columns.", file=sys.stderr)
        sys.exit(1)

    for col in numeric_cols:
        if col.lower() == "concentration" or col.lower() == "sample": continue
        med, lower, upper = median_ci(df[col], n_boot=n_boot)
        results.append((col, med, lower, upper))

    res_df = pd.DataFrame(results, columns=["Features", "Median", "CI_Lower", "CI_Upper"])
    # Round numeric columns to 4 decimals
    res_df[["Median", "CI_Lower", "CI_Upper"]] = res_df[["Median", "CI_Lower", "CI_Upper"]].round(4)

    # Save and show
    try:
        res_df.to_csv(output_csv, index=False)
    except Exception as e:
        print(f"Error: failed to write '{output_csv}': {e}", file=sys.stderr)
        sys.exit(1)

    print(res_df)

if __name__ == "__main__":
    main()
