import os
import argparse
import logging
import pickle
import tqdm

import openmm
from openmm import app, unit
from openmmtools.integrators import LangevinIntegrator

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read filename
parser = argparse.ArgumentParser(description='run equil')
parser.add_argument('outdir', type=str, help='path to files')
parser.add_argument('phase', type=str, help='apo or complex')
parser.add_argument('endstate', type=int, help='0 or 1')
parser.add_argument('length', type=float, help='in ns')
args = parser.parse_args()

# Load htf
dir_num = os.path.basename(os.path.dirname(args.outdir))
with open(os.path.join(args.outdir, f"{dir_num}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system
positions = htf.hybrid_positions

# Set parameters
_logger.info("Setting parameters...")
temperature = 300 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 4.0 * unit.femtoseconds  # We can use a 4fs timestep with HMR
platform_name = 'CUDA'

# Set number of steps
nsteps_per_iteration = 250
nequil = int(args.length * 1000)

# Load htf
integrator = LangevinIntegrator(temperature=temperature, timestep=timestep, collision_rate=collision_rate)
platform = openmm.Platform.getPlatformByName(platform_name)

# Create context
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Set context parameters
context_parameters = context.getParameters().keys()
for param in context_parameters:
    if 'old' in param:
        context.setParameter(param, 1 - args.endstate)
        _logger.info(f"Setting {param} to {1 - args.endstate}")

    elif 'new' in param:
        context.setParameter(param, args.endstate)
        _logger.info(f"Setting {param} to {args.endstate}")

# Minimize
_logger.info(f'Starting to minimize')
openmm.LocalEnergyMinimizer.minimize(context)

# Equilibrate
_logger.info(f'Starting to equilibrate for {nequil*nsteps_per_iteration*timestep}')
positions, energies = list(), list()
for _ in tqdm.tqdm(nequil):
    integrator.step(nsteps_per_iteration)
    state = context.getState(getEnergy=True, getPositions=True)
    positions.append(state.getPositions())
    energies.append(state.getPotentialEnergy())
_logger.info(f'Relax done')

# Save positions and energies
with open(os.path.join(args.outdir, f"{dir_num}_{phase}_positions.pickle"), 'wb') as f:
    pickle.dump(positions, f)

with open(os.path.join(args.outdir, f"{dir_num}_{phase}_energies.pickle"), 'wb') as f:
    pickle.dump(energies, f)


