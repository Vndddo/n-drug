#STEP 1: Import libraries
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
mol = Chem.MolFromSmiles("CCO")
print(mol)

#STEP 2: Define datasets
data = {
    'name' : ['Aspirin', 'Caffeine', 'Amphetamine', 'Ibuprofen'],
    'smiles': [
        'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)',  # Caffeine
        'CC(CC1=CC=CC=C1)N',  # Amphetamine
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'  # Ibuprofen
        ]
    }

df = pd.DataFrame(data) #turn data into DataFrame table via panda

#STEP 3: Lipinski calculation function
def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to molecule object
    mw = Descriptors.MolWt(mol)       # Molecular weight
    logp = Descriptors.MolLogP(mol)   # LogP: fat solubility
    hbd = Descriptors.NumHDonors(mol) # Number of H-bond donors
    hba = Descriptors.NumHAcceptors(mol) # Number of H-bond acceptors
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return pd.Series([mw, logp, hbd, hba, violations])

#STEP 4: Function application to data
df[['MolWt', 'LogP', 'HBD', 'HBA', 'Violations']] = df['smiles'].apply(calculate_lipinski)

#STEP 5: Drug-like compounds flag
df['Drug-like'] = df['Violations'] <= 1

#STEP 6: Print results
print(df)

#STEP 7: Saving results
df.to_csv("drug_likeness_results.csv", index=False)
