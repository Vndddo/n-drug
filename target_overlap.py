# Importing Library
import pandas as pd
import os
import re

# Configuring Input
omim_file = "OMIM-Entry-Retrieval.xlsx"
genecards_file = "GeneCards-SearchResults.xlsx"
output_file = "target_overlap_results.xlsx"

# File Configuration
def load_file(file_path, skiprows=0):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        return pd.read_csv(file_path, skiprows=skiprows)
    elif ext in [".xls", ".xlsx"]:
        return pd.read_excel(file_path, skiprows=skiprows)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

# Loading Files
omim_df = pd.read_excel(omim_file, skiprows=4)
print("OMIM Columns:", omim_df.columns.tolist())
genecards_df = pd.read_excel(genecards_file)
print("GeneCards Columns:", genecards_df.columns.tolist())

# Extracting Data From OMIM
def extract_gene_symbol(title):
    match = re.search(r';\s*(([A-Z0-9\-]+)$)', str(title))
    return match.group(1) if match else None

omim_df["Gene Symbol"] = omim_df["Title"].apply(extract_gene_symbol)

# Cleaning GeneCards Data
genecards_df["Gene Symbol"] = genecards_df["Gene Symbol"].str.strip().str.upper()

# Analysing Overlaps
overlap_df = pd.merge(genecards_df, omim_df, on="Gene Symbol", how="inner")

# Saving Results
with pd.ExcelWriter(output_file) as writer:
    omim_df.to_excel(writer, sheet_name="OMIM", index=False)
    genecards_df.to_excel(writer, sheet_name="GeneCards", index=False)
    overlap_df.to_excel(writer, sheet_name="Overlap", index=False)

print(f"Done! Results saved to: {output_file}")
print(f"Overlapping genes found: {len(overlap_df)}")
