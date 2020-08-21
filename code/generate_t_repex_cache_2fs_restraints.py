import argparse
import math
from simtk import unit
from openmmtools import testsystems, states, mcmc, multistate
import os
import tempfile
import pickle
from perses.annihilation.lambda_protocol import RelativeAlchemicalState
import mdtraj as md
import numpy as np
from simtk.unit.quantity import Quantity

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
args = parser.parse_args()

dir_num = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{dir_num}_{args.phase}.pickle"), "rb" ))

# # Add RMSD force
# from simtk.openmm import RMSDForce
# rmsd_force = RMSDForce(htf.hybrid_positions, [0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 19, 20, 21])
# htf.hybrid_system.addForce(rmsd_force)

# Multiply force constant in PeriodicTorsionForce by 100 for heavy atom non-sidechain dihedrals
atom_indices = htf.hybrid_topology.select("not name hydrogen and not sidechain")
force = htf.hybrid_system.getForce(5)
for i in range(force.getNumTorsions()):
    torsion = force.getTorsionParameters(i)
    atoms = torsion[:4]
    result = all(atom in atom_indices for atom in atoms) # Check that all atom indices are in non-sidechain heavy atom list
    if result:
        print(i, torsion)
        force.setTorsionParameters(i, torsion[0], torsion[1], torsion[2], torsion[3], torsion[4], torsion[5], torsion[6]*10)
print(htf.hybrid_system.getForce(5).getTorsionParameters(47))

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(j) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for j in range(n_replicas)]

alchemical_state = RelativeAlchemicalState.from_system(htf.hybrid_system)
alchemical_state.set_alchemical_parameters(args.state)

thermodynamic_states = [states.ThermodynamicState(system=htf.hybrid_system, temperature=T)
                        for T in temperatures]
compound_states = [states.CompoundThermodynamicState(thermodynamic_state=thermodynamic_state, composable_states=[alchemical_state])
                        for thermodynamic_state in thermodynamic_states]

# Set up sampler
move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=500)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{dir_num}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
print(reporter_file)
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=states.SamplerState(htf.hybrid_positions),
                  storage=reporter)
simulation.run()