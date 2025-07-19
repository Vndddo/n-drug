# n-drug
# Natural Product Drug Discovery Toolkit

This repository contains a suite of Python scripts designed to support **computational natural product drug discovery** through compound structure retrieval, Lipinski filtering, and target mapping.

---

## Project Overview

Natural products are a valuable source of therapeutic leads. This toolkit helps researchers streamline three essential stages in early-phase natural product screening:

1. **Structural Data Retrieval from PubChem**  
2. **Drug-likeness Screening via Lipinski's Rule of Five**  
3. **Target Prediction Overlap between OMIM and GeneCards**  

---

## Included Scripts

### `get_pubchem_data.py`
- **Purpose**: Automates retrieval of IUPAC names, molecular formulae, and SMILES strings from PubChem.
- **Input**: Excel file containing compound names.
- **Output**: Excel file with enriched compound information.
- **Usage**:
  ```bash
  python get_pubchem_data.py

### rdkit_lipinski.py
- **Purpose**: Calculates Lipinski’s drug-likeness descriptors from SMILES.
- **Input**: Excel file with SMILES column (output of get_pubchem_data.py).
- **Output**: CSV with MW, LogP, H-bond donor/acceptors, and drug-likeness flags.
- **Usage**:
  ```bash
  python rdkit_lipinski.py

### target_overlap.py
- **Purpose**: Finds overlapping gene targets between OMIM and GeneCards databases.
- **Input**: OMIM and GeneCards data files (Excel or CSV).
- **Output**: Excel file with separate sheets for OMIM, GeneCards, and overlapping targets.
- **Usage**:
  ```bash
  python target_overlap.py

## Applications
1. Natural product-based lead compound identification
2. Predictive target mapping from secondary metabolites
3. Early-phase drug-likeness filtering

## Requirements

Install dependencies via conda:
 ```bash
  conda create -n np-discovery python=3.10
  conda activate np-discovery
  conda install -c conda-forge rdkit pandas openpyxl

