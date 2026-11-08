import os
os.environ["MPLBACKEND"] = "Agg"  # Must be set BEFORE importing matplotlib in Colab subprocesses
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
import numpy as np

traj_path = "results/trajectory.dcd"
top_path = "results/equilibrated_gpr32.pdb"
plddt_path = "data/af3_plddt.csv"

# 1. Load trajectory and isolate protein atoms
traj = md.load(traj_path, top=top_path)
protein_indices = traj.topology.select("protein")
traj_prot = traj.atom_slice(protein_indices)

# Select Full-Length C-alpha (residues 1-356) vs. 7TM Helical Core C-alpha (residues 36-335)
ca_all = traj_prot.topology.select("name CA")
ca_core = traj_prot.topology.select("name CA and (resi 35 to 334)")
core_atom_indices = traj_prot.topology.select("resi 35 to 334")

# Superpose trajectory onto the stable 7TM core
traj_prot.superpose(traj_prot, frame=0, atom_indices=ca_core)

# 2. Compute Full-Length vs. 7TM Core RMSD and Radius of Gyration (in Angstroms)
rmsd_all = md.rmsd(traj_prot, traj_prot, frame=0, atom_indices=ca_all) * 10.0
rmsd_core = md.rmsd(traj_prot, traj_prot, frame=0, atom_indices=ca_core) * 10.0

rg_all = md.compute_rg(traj_prot) * 10.0
rg_core = md.compute_rg(traj_prot.atom_slice(core_atom_indices)) * 10.0
time_ps = np.arange(len(rmsd_all)) * 2.0

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 6.5), sharex=True)
ax1.plot(time_ps, rmsd_all, color="#7f7f7f", linestyle="--", alpha=0.8, lw=1.5, label="Full-Length Cα RMSD (1–356, incl. flexible tails)")
ax1.plot(time_ps, rmsd_core, color="#1f77b4", lw=2.2, label="7TM Helical Core Cα RMSD (36–335)")
ax1.set_ylabel("Cα RMSD (Å)", fontweight="bold")
ax1.set_title("Cα RMSD and radius of gyration (restrained run)", fontweight="bold")
ax1.legend(loc="upper left", frameon=True)
ax1.grid(alpha=0.3)

ax2.plot(time_ps, rg_all, color="#7f7f7f", linestyle="--", alpha=0.8, lw=1.5, label="Full-Length Rg (1–356)")
ax2.plot(time_ps, rg_core, color="#2ca02c", lw=2.2, label="7TM Helical Core Rg (36–335)")
ax2.set_xlabel("Time (ps)", fontweight="bold")
ax2.set_ylabel("Radius of Gyration (Å)", fontweight="bold")
ax2.legend(loc="upper left", frameon=True)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("results/rmsd_rg_profile.png", dpi=300)
plt.close()

# 3. Compute C-alpha RMSF (aligned on 7TM Core) and compare against AF3 pLDDT
rmsf_angstrom = md.rmsf(traj_prot, traj_prot, frame=0, atom_indices=ca_all) * 10.0
df_plddt = pd.read_csv(plddt_path)

n_res = min(len(rmsf_angstrom), len(df_plddt))
res_ids = df_plddt["residue_index"].values[:n_res]
plddt_vals = df_plddt["plddt"].values[:n_res]
rmsf_vals = rmsf_angstrom[:n_res]

# Exclude residues 1-3 to prevent the extreme terminal spike from flattening the plot scale
plot_mask = res_ids > 3
res_ids_filtered = res_ids[plot_mask]
rmsf_vals_filtered = rmsf_vals[plot_mask]
plddt_vals_filtered = plddt_vals[plot_mask]

fig, ax1 = plt.subplots(figsize=(10.5, 4.8))
ax1.axvspan(4, 35, color="#ffbb78", alpha=0.25, label="Disordered N-Terminus (4–35, residues 1–3 excluded)")
ax1.axvspan(175, 205, color="#98df8a", alpha=0.25, label="ECL2 Pocket-Gating Loop (175–205)")
ax1.axvspan(336, 356, color="#ffbb78", alpha=0.25, label="Disordered C-Terminus (336–356)")

ax1.plot(res_ids_filtered, rmsf_vals_filtered, color="#d62728", label="MD Cα RMSF (Å)", lw=1.8)
ax1.set_xlabel("Residue Number", fontweight="bold")
ax1.set_ylabel("MD Cα RMSF (Å)", color="#d62728", fontweight="bold")
ax1.tick_params(axis="y", labelcolor="#d62728")
ax1.grid(alpha=0.3)

ax2 = ax1.twinx()
ax2.plot(res_ids_filtered, plddt_vals_filtered, color="#1f77b4", alpha=0.8, linestyle="--", label="AF3 pLDDT", lw=1.5)
ax2.set_ylabel("AlphaFold 3 pLDDT Score", color="#1f77b4", fontweight="bold")
ax2.tick_params(axis="y", labelcolor="#1f77b4")
ax2.set_ylim(25, 100)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
# Move legend to lower right to avoid overlapping curves
ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower right", ncol=2, fontsize=8.5, frameon=True)

plt.title("AlphaFold 3 pLDDT vs MD Cα RMSF", fontweight="bold")
plt.tight_layout()
plt.savefig("results/plddt_vs_rmsf.png", dpi=300)
plt.close()

# 4. Save per-residue summary table
domains = [
    "N-Terminus (Disordered)" if r <= 35 else
    "C-Terminus (Disordered)" if r >= 336 else
    "ECL2 Loop" if 175 <= r <= 205 else
    "7TM Core / Short Loops"
    for r in res_ids
]

df_out = pd.DataFrame({
    "residue_index": res_ids,
    "residue_name": df_plddt["residue_name"].values[:n_res],
    "domain_region": domains,
    "af3_plddt": plddt_vals,
    "md_rmsf_angstrom": np.round(rmsf_vals, 3)
})
df_out.to_csv("results/residue_plddt_rmsf.csv", index=False)

corr = np.corrcoef(rmsf_vals, plddt_vals)[0, 1]
print(f"[INFO] Mean 7TM Core RMSD (36-335): {rmsd_core.mean():.2f} Å | Mean Full-Length RMSD (1-356): {rmsd_all.mean():.2f} Å")
print(f"[INFO] Pearson correlation between AF3 pLDDT and MD RMSF: r = {corr:.3f}")
