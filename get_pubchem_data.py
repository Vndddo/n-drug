import pandas as pd
import requests
import time

# Loading the Source file
input_file = "compounds.xlsx"  # Change to .csv if needed
df = pd.read_excel(input_file)  # or pd.read_csv for CSV
compound_names = df["Name"].dropna().tolist()

# Fetching data from PubChem
def fetch_pubchem_info(name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    try:
        # Get CID
        cid_url = f"{base_url}/compound/name/{name}/cids/JSON"
        cid_res = requests.get(cid_url, timeout=10)
        cid_res.raise_for_status()
        cid = cid_res.json()["IdentifierList"]["CID"][0]

        # Get compound properties
        prop_url = f"{base_url}/compound/cid/{cid}/property/IUPACName,MolecularFormula,CanonicalSMILES/JSON"
        prop_res = requests.get(prop_url, timeout=10)
        prop_res.raise_for_status()
        props = prop_res.json()["PropertyTable"]["Properties"][0]

        return {
            "Name": name,
            "IUPAC": props.get("IUPACName"),
            "Formula": props.get("MolecularFormula"),
            "SMILES": props.get("SMILES")
        }
    except Exception as e:
        return {"Name": name, "IUPAC": None, "Formula": None, "SMILES": None}

# Runing query
results = []
for name in compound_names:
    print(f"Processing: {name}")
    results.append(fetch_pubchem_info(name))
    time.sleep(0.2)

# Saving results to Excel
output_df = pd.DataFrame(results)
output_df.to_excel("compound_info_output.xlsx", index=False)

print("Done. Data saved to 'compound_info_output.xlsx'")
