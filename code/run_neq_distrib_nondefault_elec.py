import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from simtk import unit
from simtk import openmm
import argparse
import os
import time

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('sim_number', type=str, help='number in job name - 1')
args = parser.parse_args()

# Define lambda protocol
custom_functions = {'lambda_sterics_core':
                         lambda x: x,
                         'lambda_electrostatics_core':
                         lambda x: x,
                         'lambda_sterics_insert':
                         lambda x: 2.0 * x if x < 0.5 else 1.0,
                         'lambda_sterics_delete':
                         lambda x: 0.0 if x < 0.5 else 2.0 * (x - 0.5),
                         'lambda_electrostatics_insert':
                         lambda x: 0.0 if x < 0.5 else np.sqrt(2.0 * (x - 0.5)),
                         'lambda_electrostatics_delete':
                         lambda x: 2.0 * x if x < 0.5 else 1.0,
                         'lambda_bonds':
                         lambda x: x,
                         'lambda_angles':
                         lambda x: x,
                         'lambda_torsions':
                         lambda x: x
                         }
lambda_protocol = LambdaProtocol(functions=custom_functions)

# Define simulation parameters
nsteps_eq = 62500 # 0.25 ns
nsteps_neq= 20000 # 80 ps
neq_splitting='V R H O R V'
timestep=4.0 * unit.femtosecond
platform_name = 'CUDA'

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system
positions = htf.hybrid_positions

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(lambda_protocol, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run neq
ncycles = 10
forward_works_master, reverse_works_master = list(), list()
for cycle in range(ncycles):
    # Equilibrium (lambda = 0)
    _logger.info(f'Cycle: {cycle}, Starting to equilibrate at lambda = 0')
    initial_time = time.time()
    integrator.step(nsteps_eq)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    _logger.info(f'Cycle: {cycle}, Done equilibrating at lambda = 0, took: {elapsed_time / unit.seconds} seconds')

    # Forward (0 -> 1)
    forward_works = [integrator.get_protocol_work(dimensionless=True)]
    for fwd_step in range(nsteps_neq):
        initial_time = time.time()
        integrator.step(1)
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Cycle: {cycle}, forward NEQ step: {fwd_step}, took: {elapsed_time / unit.seconds} seconds')
        forward_works.append(integrator.get_protocol_work(dimensionless=True))
    forward_works_master.append(forward_works)

    # Equilibrium (lambda = 1)
    _logger.info(f'Cycle: {cycle}, Starting to equilibrate at lambda = 1')
    initial_time = time.time()
    integrator.step(nsteps_eq)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    _logger.info(f'Cycle: {cycle}, Done equilibrating at lambda = 1, took: {elapsed_time / unit.seconds} seconds')

    # Reverse work (1 -> 0)
    reverse_works = [integrator.get_protocol_work(dimensionless=True)]
    for rev_step in range(nsteps_neq):
        initial_time = time.time()
        integrator.step(1)
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Cycle: {cycle}, reverse NEQ step: {rev_step}, took: {elapsed_time / unit.seconds} seconds')
        reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    reverse_works_master.append(reverse_works)
        
# Save works
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward.npy"), 'wb') as f:
    np.save(f, forward_works_master)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse.npy"), 'wb') as f:
    np.save(f, reverse_works_master)