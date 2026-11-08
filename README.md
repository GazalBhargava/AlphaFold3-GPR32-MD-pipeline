# AlphaFold 3 Model and Short MD Simulation of Human GPR32

## Overview

GPR32 (UniProt O75388) is a 356-residue Class A GPCR, also known as the resolvin D1 receptor. To my knowledge, no experimental structure of human GPR32 has been deposited in the RCSB PDB (checked as of November 2026). This project models it with AlphaFold 3 and runs a short molecular dynamics (MD) simulation to see which parts of the model move. It is an exploratory exercise, not a validated structure.

## Results

**AlphaFold 3 model.** The model (AlphaFold Server, single chain) has pTM 0.78. Mean pLDDT is 86.3 over residues 36–335, which I took as the seven-transmembrane (7TM) core, and 78.1 over the full chain. The tails are low confidence: mean pLDDT 28.5 for the N-terminus (1–35) and 43.0 for the C-terminus (336–356). Inside the core, the lowest-confidence stretches (pLDDT below 70) are residues 196–207, 249–255, 317–324 and 331–335. In extracellular loop 2 (ECL2, 175–205), the first half (175–195) is well predicted (mean pLDDT 83.6) and the tip (196–205) is not (mean 55.3).

**MD simulation.** During the 500 ps production run, the 219 residues with pLDDT of 85 or higher were restrained and stayed within 0.72 Å RMSF (mean 0.47 Å). That is set by the restraints and does not measure stability. The unrestrained regions moved much more: the ECL2 tip (196–205) had a mean RMSF of 2.7 Å, the C-terminus 4.7 Å and the N-terminus (residues 4–35) 5.7 Å. Residues 1–3 are left out of these averages because their RMSF (51–58 Å) is far above that of residue 4 (10 Å), which cannot be right for bonded neighbours, so I treat them as unreliable. The temperature averaged 310.4 K (SD 3.4 K), and the potential energy fell by about 700 kJ/mol between the first and last 50 ps, so the system was still relaxing.

**pLDDT and fluctuation.** Among the 137 unrestrained residues, lower pLDDT goes with higher RMSF (Spearman ρ = −0.94). Across all 356 residues ρ = −0.87, but that figure is partly built in, because the restraints were assigned from pLDDT.

## Simulation setup

| Setting | Value |
| --- | --- |
| Starting model | AlphaFold 3 mmCIF, protonated at pH 7.4 with PDBFixer |
| Force field and solvent | Amber14 with GBn2 implicit solvent. No explicit water, ions or membrane |
| Restraints | Harmonic position restraints (k = 500 kJ/mol/nm²) on the Cα atoms of residues with pLDDT of 85 or higher, as a crude stand-in for the membrane |
| Integrator | Langevin middle integrator, 310 K, 2 fs time step, bonds to hydrogen constrained |
| Protocol | 1,000 minimisation iterations, 20 ps equilibration, 500 ps production (one trajectory, a frame every 2 ps) |
| Hardware | OpenMM on a Google Colab T4 GPU (High-RAM runtime) |

Per-residue Cα RMSF is measured against the first production frame after superposing on residues 36–335.

## Limitations

- One 500 ps run in implicit solvent without a membrane is too short and too simple to sample loop conformations or give error bars.
- The restrained helices cannot show that the 7TM core is intrinsically stable.
- RMSF against the starting frame mixes drift with fluctuation.
- There is no ligand, G protein or experimental comparison.

## Next steps

- Embed the receptor in an explicit POPC bilayer with water and ions, release the restraints and run longer replicas.
- Sample the ECL2 tip more thoroughly to generate open and closed pocket states, then dock resolvin D1 and lipoxin A4.

## Running the workflow

The notebook `af3_gpr32_md_pipeline.ipynb` runs the whole workflow on Google Colab with a GPU: it extracts pLDDT from the AlphaFold Server model, protonates it, runs the simulation and produces the plots and the per-residue table. It needs OpenMM, PDBFixer, MDTraj, Biopython, pandas, NumPy and Matplotlib.
