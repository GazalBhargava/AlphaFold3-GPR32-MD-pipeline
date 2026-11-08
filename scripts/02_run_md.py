import os, sys
import pandas as pd
from openmm.app import *
from openmm import *
from openmm.unit import *

os.makedirs("results", exist_ok=True)
pdb = PDBFile("data/gpr32_clean.pdb")

# 1. Implicit-solvent GPCR system (Amber14 + GBn2)
forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
system = forcefield.createSystem(
    pdb.topology, nonbondedMethod=CutoffNonPeriodic,
    nonbondedCutoff=2.0*nanometer, constraints=HBonds
)

# 2. Anchor 7TM helices (AF3 pLDDT >= 85) to mimic membrane; keep loops/termini free
df_plddt = pd.read_csv("data/af3_plddt.csv")
tm_residues = set(df_plddt[df_plddt["plddt"] >= 85.0]["residue_index"].astype(int))

restraint = CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
restraint.addGlobalParameter("k", 500.0 * kilojoules_per_mole / nanometer**2)
for p in ("x0", "y0", "z0"):
    restraint.addPerParticleParameter(p)

for atom in pdb.topology.atoms():
    if atom.name == "CA" and int(atom.residue.id) in tm_residues:
        restraint.addParticle(atom.index, pdb.positions[atom.index].value_in_unit(nanometer))
system.addForce(restraint)

# 3. Minimize & equilibrate (OpenMM auto-selects CUDA GPU)
integrator = LangevinMiddleIntegrator(310*kelvin, 1.0/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(maxIterations=1000)
simulation.context.setVelocitiesToTemperature(310*kelvin)
simulation.step(10000)

with open("results/equilibrated_gpr32.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

# 4. Production MD (250,000 steps = 500 ps)
simulation.reporters.extend([
    DCDReporter("results/trajectory.dcd", 1000),
    StateDataReporter("results/md_log.csv", 1000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True),
    StateDataReporter(sys.stdout, 25000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True)
])
simulation.step(250000)
