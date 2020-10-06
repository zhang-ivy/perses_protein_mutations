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
from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping


# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('direction', type=str, help="forward or reverse")
args = parser.parse_args()

i = os.path.basename(os.path.dirname(args.dir))

with open(os.path.join(args.dir, f"{i}_vacuum.pickle" ), "rb") as f:
    htf_at = pickle.load(f)

old_system = htf_at._topology_proposal.old_system
new_system = htf_at._topology_proposal.new_system

nb_force = old_system.getForce(3)
for i in range(nb_force.getNumExceptions()):
    p1, p2, chargeProd, sigma, epsilon = nb_force.getExceptionParameters(i)
    if p1 in [10, 11, 12, 13] or p2 in [10, 11, 12, 13]:
        nb_force.setExceptionParameters(i, p1, p2, 0, sigma, 0)

nb_force = new_system.getForce(3)
for i in range(nb_force.getNumExceptions()):
    p1, p2, chargeProd, sigma, epsilon = nb_force.getExceptionParameters(i)
    if p1 in [10, 13, 14, 15, 16, 17, 18, 19] or p2 in [10, 13, 14, 15, 16, 17, 18, 19]:
        nb_force.setExceptionParameters(i, p1, p2, 0, sigma, 0)

# Build new htf
htf_at._topology_proposal._old_system = old_system
htf_at._topology_proposal._new_system = new_system

from perses.annihilation.relative import HybridTopologyFactory
htf = HybridTopologyFactory(topology_proposal=htf_at._topology_proposal,
                     current_positions=htf_at.old_positions(htf_at.hybrid_positions),
                     new_positions=htf_at.new_positions(htf_at.hybrid_positions),
                     use_dispersion_correction=False,
                     functions=None,
                     softcore_alpha=None,
                     bond_softening_constant=1.0,
                     angle_softening_constant=1.0,
                     soften_only_new=False,
                     neglected_new_angle_terms=[],
                     neglected_old_angle_terms=[],
                     softcore_LJ_v2=True,
                     softcore_electrostatics=True,
                     softcore_LJ_v2_alpha=0.85,
                     softcore_electrostatics_alpha=0.3,
                     softcore_sigma_Q=1.0,
                     interpolate_old_and_new_14s=False,
                     omitted_terms=None,
                     flatten_torsions=True)

pickle.dump(htf, open(os.path.join(args.dir, f"{i}_vacuum_{args.name.lower()}.pickle"), "wb" ))

# Render atom map
atom_map_filename = f'{args.dir}/atom_map_{args.name.lower()}.png'
render_protein_residue_atom_mapping(htf._topology_proposal, atom_map_filename)

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

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
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=states.SamplerState(htf.hybrid_positions),
                  storage=reporter)
simulation.run()