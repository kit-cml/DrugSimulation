import pandas as pd
import numpy as np
import glob
import os

RESULT_ROOT = "results"

# Find all dose folders
dose_folders = [ 
    f.path for f in os.scandir(RESULT_ROOT) if f.is_dir()
]

if len(dose_folders) == 0:
    raise RuntimeError("No dose folders found inside results/")

output_rows = []

for folder in sorted(dose_folders):
    dose_name = os.path.basename(folder)

    csv_files = glob.glob(os.path.join(folder, "*features*.csv"))

    if len(csv_files) == 0:
        print("Skipping (no CSV files):", folder)
        continue

    print("Processing dose:", dose_name)

    dfs = []
    for f in csv_files:
        df = pd.read_csv(f)

        # Drop the sample column if it exists
        if "sample" in df.columns:
            df = df.drop(columns=["sample"])

        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)

    # Replace weird strings with real NaN
    combined = combined.replace(['-nan', 'nan', 'NaN'], np.nan)

    # Compute column means
    mean_vals = combined.mean(numeric_only=True)

    # Store result
    row = {"dose": dose_name}
    row.update(mean_vals.to_dict())
    output_rows.append(row)

# -------------------------------------------------------------------------
# Convert list of dict → DataFrame
# -------------------------------------------------------------------------
result_df = pd.DataFrame(output_rows)

# Convert dose to numeric and sort
result_df["dose"] = result_df["dose"].astype(float)
result_df = result_df.sort_values(by="dose")

# Round numeric values
numeric_cols = result_df.select_dtypes(include=[float, np.number]).columns
result_df[numeric_cols] = result_df[numeric_cols].round(4)

# Save normal table
output_csv = os.path.join(RESULT_ROOT, "dose_averages.csv")
result_df.to_csv(output_csv, index=False)
print("\nAverages written to dose_averages.csv")

# -------------------------------------------------------------------------
# Create TRANSPOSED RESULT (features = rows, doses = columns)
# -------------------------------------------------------------------------
transposed = result_df.set_index("dose").transpose()

# Label first column
transposed.index.name = "Features---Concentrations (nM)"

# Save transposed CSV
output_csv_t = os.path.join(RESULT_ROOT, "dose_averages_transposed.csv")
transposed.to_csv(output_csv_t)

# Save transposed Markdown (compatible with older pandas)
output_md_t = os.path.join(RESULT_ROOT, "dose_averages_transposed.md")

with open(output_md_t, "w") as f:
    f.write(transposed.to_markdown())

print("Transposed averages written to dose_averages_transposed.csv and .md")
