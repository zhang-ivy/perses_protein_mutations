import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from simtk import unit
from simtk import openmm
import argparse
import os
import time
import mdtraj as md

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('sim_number', type=str, help='number in job name - 1')
args = parser.parse_args()

# Define lambda protocol
x = 'lambda'
custom_functions = {
                     'lambda_sterics_core': x,
                     'lambda_electrostatics_core': x,
                     'lambda_sterics_insert': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                     'lambda_sterics_delete': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                     'lambda_electrostatics_insert': f"select(step({x} - 0.5), sqrt(2.0 * ({x} - 0.5)), 0.0)",
                     'lambda_electrostatics_delete': f"select(step({x} - 0.5), 1.0, sqrt(2.0 * {x}))",
                     'lambda_bonds': x,
                     'lambda_angles': x,
                     'lambda_torsions': x}

# Define simulation parameters
nsteps_eq = 62500 # 0.25 ns
nsteps_neq = 20000 # 80 ps
neq_splitting ='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system
positions = htf.hybrid_positions

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(custom_functions, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep)

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
forward_traj_old, forward_traj_new, reverse_traj_old, reverse_traj_new = list(), list(), list(), list()
for cycle in range(ncycles):
    # Equilibrium (lambda = 0)
    _logger.info(f'Cycle: {cycle}, Starting to equilibrate at lambda = 0')
    initial_time = time.time()
    integrator.step(nsteps_eq)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    _logger.info(f'Cycle: {cycle}, Done equilibrating at lambda = 0, took: {elapsed_time / unit.seconds} seconds')

    # Forward (0 -> 1)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    forward_traj_old.append(old_pos)
    forward_works = [integrator.get_protocol_work(dimensionless=True)]
    for fwd_step in range(nsteps_neq):
        initial_time = time.time()
        integrator.step(1)
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Cycle: {cycle}, forward NEQ step: {fwd_step}, took: {elapsed_time / unit.seconds} seconds')
        forward_works.append(integrator.get_protocol_work(dimensionless=True))
    forward_works_master.append(forward_works)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    new_pos = np.asarray(htf.new_positions(pos))
    forward_traj_new.append(new_pos)

    # Equilibrium (lambda = 1)
    _logger.info(f'Cycle: {cycle}, Starting to equilibrate at lambda = 1')
    initial_time = time.time()
    integrator.step(nsteps_eq)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    _logger.info(f'Cycle: {cycle}, Done equilibrating at lambda = 1, took: {elapsed_time / unit.seconds} seconds')

    # Reverse work (1 -> 0)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    new_pos = np.asarray(htf.new_positions(pos))
    reverse_traj_new.append(new_pos)
    reverse_works = [integrator.get_protocol_work(dimensionless=True)]
    for rev_step in range(nsteps_neq):
        initial_time = time.time()
        integrator.step(1)
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Cycle: {cycle}, reverse NEQ step: {rev_step}, took: {elapsed_time / unit.seconds} seconds')
        reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    reverse_works_master.append(reverse_works)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    reverse_traj_old.append(old_pos)
        
        
# Save works
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward.npy"), 'wb') as f:
    np.save(f, forward_works_master)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse.npy"), 'wb') as f:
    np.save(f, reverse_works_master)

# Save positions
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_old.npy"), 'wb') as f:
    np.save(f, forward_traj_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_new.npy"), 'wb') as f:
    np.save(f, forward_traj_new)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_old.npy"), 'wb') as f:
    np.save(f, reverse_traj_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_new.npy"), 'wb') as f:
    np.save(f, reverse_traj_new)