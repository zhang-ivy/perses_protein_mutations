
from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState
from openmmtools import states
from simtk import openmm, unit
from openmmtools.mcmc import LangevinSplittingDynamicsMove, GHMCMove
from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import numpy as np
import pickle
import argparse
import os

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA, corresponding to the endstate')
parser.add_argument('endstate', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
args = parser.parse_args()

# Read htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.endstate}.pickle"), "rb") as f:
    htf = pickle.load(f)

# Alchemify the hybrid system
# atoms_to_alchemify = list(htf._atom_classes['unique_new_atoms']) + list(htf._atom_classes['unique_old_atoms'])
if args.endstate == 0:
    atoms_to_alchemify = list(htf._atom_classes['unique_old_atoms'])
elif args.endstate == 1:
    atoms_to_alchemify = list(htf._atom_classes['unique_new_atoms'])

alch_factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
alchemical_region = AlchemicalRegion(alchemical_atoms=list(atoms_to_alchemify), alchemical_torsions=True)
alchemical_system = alch_factory.create_alchemical_system(htf.hybrid_system, alchemical_region)

# Initialize compound thermodynamic states at different temperatures and alchemical states.
protocol = {'temperature': [300]*unit.kelvin*5,
            'lambda_electrostatics': [1.0, 0.1, 0.01, 0.001, 0],
            'lambda_sterics': [1.0, 0.1, 0.01, 0.001, 0],
           'lambda_torsions': [1.0, 0.1, 0.01, 0.001, 0]}

alchemical_state = AlchemicalState.from_system(alchemical_system)
compound_states = states.create_thermodynamic_state_protocol(alchemical_system, 
                                                             protocol=protocol, 
                                                             composable_states=[alchemical_state])

# Set up the sampler
n_steps = 500 # 1 ps
n_iterations = args.length*1000 

# # Propagate the replicas with Langevin dynamics.
# langevin_move = LangevinSplittingDynamicsMove(timestep=2.0*unit.femtosecond, n_steps=n_steps)
# simulation = ReplicaExchangeSampler(mcmc_moves=langevin_move, number_of_iterations=n_iterations)

# Propagate the replicas with GHMC move.
ghmc_move = GHMCMove(timestep=4.0*unit.femtosecond, n_steps=n_steps)
simulation = ReplicaExchangeSampler(mcmc_moves=ghmc_move, number_of_iterations=n_iterations)

# Create sampler state
sampler_state = states.SamplerState(htf.hybrid_positions)
sampler_state.box_vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()

# Run simulation
i = os.path.basename(os.path.dirname(args.dir))
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=sampler_state,
                  storage=reporter)
simulation.run()
