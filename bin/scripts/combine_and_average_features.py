#!/usr/bin/env python3
import os
import argparse
import pandas as pd


def combine_and_average(
    base_dir: str,
    pattern_keyword: str = "features",               # only pick files whose names contain this
    output_suffix: str = "-features-all.csv",        # suffix for combined files (hyphen style)
    sample_avg_filename: str = "sample-averages-across-concentrations.csv",
):
    # 1) Combine all feature files within each concentration folder
    for conc in os.listdir(base_dir):
        conc_path = os.path.join(base_dir, conc)
        if not os.path.isdir(conc_path):
            continue

        prefix = None
        files = []

        # Pick only CSV files whose names contain "features"
        for fname in os.listdir(conc_path):
            if pattern_keyword in fname and fname.endswith(".csv"):
                files.append(fname)
                if prefix is None:
                    # Example: quinidine_6474.00_features_core0.csv
                    # 1) Split at "features" and trim trailing underscore
                    raw_prefix = fname.split(pattern_keyword)[0].rstrip("_")
                    # 2) Replace the FIRST underscore (between drug and conc) by hyphen:
                    #    quinidine_6474.00 -> quinidine-6474.00
                    prefix = raw_prefix.replace("_", "-", 1)

        if not files or prefix is None:
            continue

        files.sort()
        dfs = [pd.read_csv(os.path.join(conc_path, f)) for f in files]
        combined = pd.concat(dfs, ignore_index=True)

        # Drop exact duplicate rows
        combined.drop_duplicates(inplace=True)

        # Sort combined file (by Sample if present, otherwise by index)
        if "Sample" in combined.columns:
            combined.sort_values("Sample", inplace=True)
        else:
            combined.sort_index(inplace=True)

        combined_out = os.path.join(conc_path, f"{prefix}{output_suffix}")
        combined.to_csv(combined_out, index=False)
        print(f"Written combined file: {combined_out}")

    # 2) Collect all combined files, tagged with their concentration
    #    Only keep directory names that can be parsed as float (concentrations)
    conc_dirs = []
    for d in os.listdir(base_dir):
        dpath = os.path.join(base_dir, d)
        if not os.path.isdir(dpath):
            continue
        try:
            float(d)
            conc_dirs.append(d)
        except ValueError:
            continue

    conc_dirs.sort(key=lambda x: float(x))

    all_rows = []
    for conc in conc_dirs:
        conc_path = os.path.join(base_dir, conc)
        combined_files = [
            f for f in os.listdir(conc_path)
            if f.endswith(output_suffix)
        ]
        if not combined_files:
            continue

        fpath = os.path.join(conc_path, combined_files[0])
        df = pd.read_csv(fpath)

        # Attach concentration info
        df["Concentration"] = float(conc)
        all_rows.append(df)

    # 3) Compute per‑sample averages across concentrations
    if all_rows:
        all_df = pd.concat(all_rows, ignore_index=True)

        # Ensure Sample exists
        if "Sample" not in all_df.columns:
            raise ValueError("Expected a 'Sample' column in the data.")

        # Numeric columns excluding Sample
        numeric_cols = [
            c for c in all_df.select_dtypes(include="number").columns
            if c != "Sample"
        ]

        per_sample_avg = (
            all_df
            .groupby("Sample")[numeric_cols]
            .mean()
        )
        per_sample_avg = per_sample_avg.reset_index()
        per_sample_avg.sort_values("Sample", inplace=True)

        sample_avg_out = os.path.join(base_dir, sample_avg_filename)
        per_sample_avg.to_csv(sample_avg_out, index=False)
        print(f"Written per-sample averages across concentrations: {sample_avg_out}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine feature CSVs per concentration, remove duplicates, "
            "and compute per-sample averages across all concentrations."
        )
    )
    parser.add_argument(
        "--base_dir",
        type=str,
        required=True,
        help="Base directory containing concentration folders (e.g. 'results' or 'results/').",
    )
    parser.add_argument(
        "--pattern_keyword",
        type=str,
        default="features",
        help="Substring that must appear in feature file names.",
    )
    parser.add_argument(
        "--output_suffix",
        type=str,
        default="-features-all.csv",
        help="Suffix for combined per‑concentration files.",
    )
    parser.add_argument(
        "--sample_avg_filename",
        type=str,
        default="sample-averages-across-concentrations.csv",
        help="Output CSV name for per-sample averages.",
    )

    args = parser.parse_args()

    combine_and_average(
        base_dir=args.base_dir,
        pattern_keyword=args.pattern_keyword,
        output_suffix=args.output_suffix,
        sample_avg_filename=args.sample_avg_filename,
    )


if __name__ == "__main__":
    main()
