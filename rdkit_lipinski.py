# Importing Libraries
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Loading Source File
input_file = "compound_info_output.xlsx"
df = pd.read_excel(input_file)

# Lipinski Calculation Function
def calculate_lipinski(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        return pd.Series([mw, logp, hbd, hba, violations])
    except Exception as e:
        print(f"Error processing SMILES: {smiles} | Error: {e}")
        return pd.Series([None, None, None, None, None])

# Apply function to SMILES column
df[['MolWt', 'LogP', 'HBD', 'HBA', 'Violations']] = df['SMILES'].apply(calculate_lipinski)

# Flag drug-likeness
df['Drug-like'] = df['Violations'] <= 1

# Save results to Excel
output_file = "drug_likeness_results.xlsx"
df.to_excel(output_file, index=False)

print(f"Analysis completed. Results saved to {output_file}")
