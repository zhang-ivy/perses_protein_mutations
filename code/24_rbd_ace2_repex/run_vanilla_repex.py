from simtk import unit, openmm

# # Load test system
# from simtk.openmm import XmlSerializer
# with open('debug_repex/16_apo_0_hybrid_system.xml', 'r') as infile:
#     alchemical_system = XmlSerializer.deserialize(infile.read())
    
# Read initial coordinates
#from simtk.openmm.app import PDBFile
#pdbfile = PDBFile('initial.pdb')
#positions = pdbfile.positions

import pickle
# with open('debug_repex/16_apo_0_hybrid_positions.pickle', 'rb') as infile:
#     positions = pickle.load(infile)


import argparse
# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
args = parser.parse_args()

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)

# Initialize compound thermodynamic states at different temperatures and alchemical states.
from openmmtools import alchemy, states
import numpy as np
n_replicas = 40
MAX_LAMBDA = 0.01
protocol = {
    'temperature': np.linspace(300, 300, n_replicas) * unit.kelvin,
    'lambda_torsions': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_angles': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_bonds': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_electrostatics_core': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_sterics_core': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_sterics_insert': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_sterics_delete': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_electrostatics_insert': np.linspace(0.0, MAX_LAMBDA, n_replicas),
    'lambda_electrostatics_delete': np.linspace(0.0, MAX_LAMBDA, n_replicas),
}
# Make sure AlchemicalState knows about some alchemical parameters
for name in ['lambda_sterics_core', 'lambda_electrostatics_core', 'lambda_sterics_delete', 'lambda_sterics_insert', 'lambda_electrostatics_delete', 'lambda_electrostatics_insert']:
    # Add the parameter
    setattr(alchemy.AlchemicalState, name, alchemy.AlchemicalState._LambdaParameter(name))
alchemical_state = alchemy.AlchemicalState.from_system(alchemical_system)
thermodynamic_states = states.create_thermodynamic_state_protocol(alchemical_system, protocol=protocol, composable_states=[alchemical_state])

# Initialize replica initial configurations.
sampler_states = [states.SamplerState(positions=positions, box_vectors=alchemical_system.getDefaultPeriodicBoxVectors()) for _ in thermodynamic_states]

# Propagate the replicas with Langevin dynamics.
from openmmtools import mcmc
n_steps = 250
#langevin_move = mcmc.LangevinSplittingDynamicsMove(timestep=2.0*unit.femtosecond, n_steps=n_steps)
langevin_move = mcmc.LangevinDynamicsMove(timestep=2.0*unit.femtosecond, n_steps=n_steps)

# Enable logging
import logging, sys
logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

# Run the parallel tempering simulation.
from openmmtools import multistate
n_iterations = 500
parallel_tempering = multistate.ReplicaExchangeSampler(mcmc_moves=langevin_move, number_of_iterations=n_iterations)
import tempfile
storage_path = f'{i}_{phase}.nc'
reporter = multistate.MultiStateReporter(storage_path, checkpoint_interval=n_iterations)
parallel_tempering.create(thermodynamic_states=thermodynamic_states, sampler_states=sampler_states, storage=reporter)
parallel_tempering.run(n_iterations=n_iterations)