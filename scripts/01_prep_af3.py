# CELL 3: Write Step 1 Script (pLDDT Extraction & Protonation)
import os
import pandas as pd
from Bio.PDB import MMCIFParser
from pdbfixer import PDBFixer
from openmm.app import PDBFile

os.makedirs("data", exist_ok=True)
cif_path = "data/af3_gpr32_model.cif"
pdb_out = "data/gpr32_clean.pdb"
plddt_out = "data/af3_plddt.csv"

# 1. Extract per-residue pLDDT (stored in the B-factor field of AF3 CIFs)
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("AF3_GPR32", cif_path)

records = []
for model in structure:
    for chain in model:
        for res in chain:
            if "CA" in res:
                records.append({
                    "residue_index": res.get_id()[1],
                    "residue_name": res.get_resname(),
                    "plddt": res["CA"].get_bfactor()
                })

df_plddt = pd.DataFrame(records)
df_plddt.to_csv(plddt_out, index=False)

core_plddt = df_plddt[(df_plddt["residue_index"] >= 36) & (df_plddt["residue_index"] <= 335)]["plddt"].mean()
print(f"[INFO] Extracted pLDDT for {len(df_plddt)} residues of GPR32.")
print(f"[INFO] Mean pLDDT (Full-Length 1-356): {df_plddt['plddt'].mean():.2f} | Mean pLDDT (7TM Core 36-335): {core_plddt:.2f}")

# 2. Clean structure and add hydrogens at physiological pH (7.4) for MD
fixer = PDBFixer(filename=cif_path)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)

with open(pdb_out, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print(f"[INFO] Saved protonated GPR32 structure to {pdb_out}")
