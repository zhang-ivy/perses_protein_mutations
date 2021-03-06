import logging
import pickle
import numpy as np
from openmmtools.integrators import LangevinIntegrator
from openmmtools.states import ThermodynamicState, CompoundThermodynamicState
from perses.annihilation.lambda_protocol import RelativeAlchemicalState, LambdaProtocol
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState
from simtk import openmm, unit
import argparse
import os
import time
import mdtraj as md
from tqdm import tqdm
from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('file', type=str, help='path to input file')
args = parser.parse_args()

# Define simulation parameters
temperature = 300 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
nsteps = 100 
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'

_logger.info("Reading in htf")
htf = pickle.load(open(args.file, "rb" ))

system = htf.hybrid_system
positions = htf.hybrid_positions
topology = htf.hybrid_topology

# # Create CompoundThermodynamicState at the appropriate endstate -- r-htf
atoms_to_alchemify = list(htf._atom_classes['unique_new_atoms']) + list(htf._atom_classes['unique_old_atoms'])
alch_factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
alchemical_region = AlchemicalRegion(alchemical_atoms=list(atoms_to_alchemify), alchemical_torsions=True, annihilate_sterics=True, annihilate_electrostatics=True)
alchemical_system = alch_factory.create_alchemical_system(system, alchemical_region)
alchemical_state = AlchemicalState.from_system(alchemical_system)
thermodynamic_state = CompoundThermodynamicState(ThermodynamicState(alchemical_system, temperature=temperature), composable_states=[alchemical_state])
thermodynamic_state.set_alchemical_variable('lambda', 1.0) # this should be 1, not the same as the endstate lambda

# Set up integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = thermodynamic_state.create_context(integrator, platform=platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run equilibration
final_pos = np.empty(shape=(101, htf.hybrid_topology.n_atoms, 3))
pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
i = 0
final_pos[i] = pos * unit.nanometers
for step in tqdm(range(nsteps)):
    i += 1
    initial_time = time.time()
    integrator.step(1)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    final_pos[i] = pos *unit.nanometers
    elapsed_time = (time.time() - initial_time) * unit.seconds
    _logger.info(f'Step: {step} took {elapsed_time} seconds')

# Save traj
with open(os.path.join(f"rbd_ace2_hybrid_pos.npy"), 'wb') as f:
    np.save(f, final_pos)
