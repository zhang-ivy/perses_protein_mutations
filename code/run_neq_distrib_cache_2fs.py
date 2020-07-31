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
parser.add_argument('old_name', type=str, help='three letter code of old amino acid')
parser.add_argument('new_name', type=str, help='three letter code of new amino acid')
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
nsteps_eq = 10
nsteps_neq = 40000 # 80 ps
neq_splitting='V R H O R V'
timestep = 2.0 * unit.femtosecond
platform_name = 'CUDA'

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system

# Read in ser cache
with open(os.path.join(args.dir, f"{args.old_name}_pos_hybrid.npy"), 'rb') as f:
    ser_pos_hybrid = np.load(f)

# Read in indices of uncorrelated ser snapshots
with open(os.path.join(args.dir, f"{args.old_name}_indices.npy"), 'rb') as f:
    ser_indices = np.load(f)

# Get equilbrium snapshot of ser
positions = ser_pos_hybrid[ser_indices[int(args.sim_number)]]

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(DEFAULT_ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep)

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

# Run eq forward (0 -> 1)
initial_time = time.time()
integrator.step(nsteps_eq)
elapsed_time = (time.time() - initial_time) * unit.seconds
_logger.info(f'Equilibrating at lambda = 0, took: {elapsed_time / unit.seconds} seconds')

# Run neq forward (0 -> 1)
forward_works_master = list()
forward_neq_old, forward_neq_new = list(), list()
forward_works = [integrator.get_protocol_work(dimensionless=True)]
for fwd_step in range(nsteps_neq):
    initial_time = time.time()
    integrator.step(1)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    forward_works.append(integrator.get_protocol_work(dimensionless=True))
    if fwd_step % 750 == 0:
        _logger.info(f'forward NEQ step: {fwd_step}, took: {elapsed_time / unit.seconds} seconds')
        pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
        old_pos = np.asarray(htf.old_positions(pos))
        new_pos = np.asarray(htf.new_positions(pos))
        forward_neq_old.append(old_pos)
        forward_neq_new.append(new_pos)
forward_works_master.append(forward_works)

# Read in ala cache
with open(os.path.join(args.dir, f"{args.new_name}_pos_hybrid.npy"), 'rb') as f:
    ala_pos_hybrid = np.load(f)

# Read in indices of uncorrelated ala snapshots
with open(os.path.join(args.dir, f"{args.new_name}_indices.npy"), 'rb') as f:
    ala_indices = np.load(f)

# Get equilbrium snapshot of ala
positions = ala_pos_hybrid[ala_indices[int(args.sim_number)]]
context.setPositions(positions)

# Run eq reverse (1 -> 0)
initial_time = time.time()
integrator.step(nsteps_eq)
elapsed_time = (time.time() - initial_time) * unit.seconds
_logger.info(f'Equilibrating at lambda = 1, took: {elapsed_time / unit.seconds} seconds')

# Run neq reverse (1 -> 0)
reverse_works_master = list()
reverse_neq_old, reverse_neq_new = list(), list()
reverse_works = [integrator.get_protocol_work(dimensionless=True)]
for rev_step in range(nsteps_neq):
    initial_time = time.time()
    integrator.step(1)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    if rev_step % 750 == 0:
        _logger.info(f'reverse NEQ step: {rev_step}, took: {elapsed_time / unit.seconds} seconds')
        pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
        old_pos = np.asarray(htf.old_positions(pos))
        new_pos = np.asarray(htf.new_positions(pos))
        reverse_neq_old.append(old_pos)
        reverse_neq_new.append(new_pos)
reverse_works_master.append(reverse_works)


# Save works
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward.npy"), 'wb') as f:
    np.save(f, forward_works_master)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse.npy"), 'wb') as f:
    np.save(f, reverse_works_master)

# Save trajs
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_old.npy"), 'wb') as f:
    np.save(f, forward_neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_new.npy"), 'wb') as f:
    np.save(f, forward_neq_new)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_old.npy"), 'wb') as f:
    np.save(f, reverse_neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_new.npy"), 'wb') as f:
    np.save(f, reverse_neq_new)

