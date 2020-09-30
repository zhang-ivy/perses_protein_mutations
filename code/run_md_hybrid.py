import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from openmmtools.states import ThermodynamicState, CompoundThermodynamicState
from perses.annihilation.lambda_protocol import RelativeAlchemicalState, LambdaProtocol
from simtk import openmm, unit
import argparse
import os
import time
from simtk.openmm.app import PDBFile
import mdtraj as md
from tqdm import tqdm

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('endstate', type=int, help="0 or 1")
args = parser.parse_args()

# Define lambda functions
x = 'lambda'
DEFAULT_ALCHEMICAL_FUNCTIONS = {
                             'lambda_sterics_core': x,
                             'lambda_electrostatics_core': x,
                             'lambda_sterics_insert': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_sterics_delete': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_insert': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_delete': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_bonds': x,
                             'lambda_angles': x,
                             'lambda_torsions': x}

# Define simulation parameters
temperature = 300 * unit.kelvin
nsteps_eq = 100000000 # 200 ns
nsteps_neq = 20000 # 80 ps
neq_splitting ='V R H O R V'
timestep = 2.0 * unit.femtosecond
platform_name = 'CUDA'

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system
positions = htf.hybrid_positions

# Create CompoundThermodynamicState at the appropriate endstate
lambda_alchemical_state = RelativeAlchemicalState.from_system(system)
lambda_protocol = LambdaProtocol(functions = 'default')
lambda_alchemical_state.set_alchemical_parameters(args.endstate, lambda_protocol)
thermodynamic_state = CompoundThermodynamicState(ThermodynamicState(system, temperature=temperature), composable_states=[lambda_alchemical_state])

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(DEFAULT_ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep, temperature=temperature)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = thermodynamic_state.create_context(integrator, platform=platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run equilibration
final_pos = np.empty(shape=(40001, htf.hybrid_topology.n_atoms, 3))
pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
i = 0
final_pos[i] = pos * unit.nanometers
for step in tqdm(range(nsteps_eq)):
    initial_time = time.time()
    integrator.step(1)
    if step % 2500 == 0:
        pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
        final_pos[i] = pos *unit.nanometers
        i += 1
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Step: {step} took {elapsed_time} seconds')

# Save traj
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}_equil_hybrid_{args.endstate}_200ns.npy"), 'wb') as f:
    np.save(f, final_pos)
