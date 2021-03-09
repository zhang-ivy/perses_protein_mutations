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
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run vanilla md')
parser.add_argument('file', type=str, help='path to input file')
parser.add_argument('--is_old', action="store_true", help='indicates we are using the new system')
args = parser.parse_args()

# Define simulation parameters
temperature = 300 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
nsteps = 100
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'

_logger.info("Reading in htf")
htf = pickle.load(open(args.file, "rb" ))

if args.is_old:
	_logger.info("Setting up old")
	system = htf._topology_proposal.old_system
	positions = htf.old_positions(htf.hybrid_positions)
	topology = htf._topology_proposal.old_topology
else:
	_logger.info("Setting up new")
	system = htf._topology_proposal.new_system
	positions = htf.new_positions(htf.hybrid_positions)
	topology = htf._topology_proposal.new_topology

# Set up integrator
_logger.info("Making integrator")
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Set up context
_logger.info("Making context")
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')

context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Minimize
_logger.info("Minimizing")
openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)

# Run equilibration
_logger.info("Equilibrating")
final_pos = np.empty(shape=(101, topology.getNumAtoms(), 3))
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
name = 'old' if args.is_old else 'new'
_logger.info(f"Saving {name}")
with open(os.path.join(f"rbd_ace2_{name}_pos.npy"), 'wb') as f:
    np.save(f, final_pos)

