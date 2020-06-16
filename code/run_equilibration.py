import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from simtk import unit
from simtk import openmm
import argparse
import os
import time
from simtk.openmm.app import PDBFile
import mdtraj as md

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
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
nsteps_eq = 625000 # 2.5 ns
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

# Run equilibration
forward_works_master, reverse_works_master = list(), list()
forward_traj_old, forward_traj_new, reverse_traj_old, reverse_traj_new = list(), list(), list(), list()

_logger.info(f'Starting to equilibrate first end state')
initial_time = time.time()
pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
old_pos = np.asarray(htf.old_positions(pos))
traj_old.append(old_pos)
for step in range(nsteps_eq):
    if step % 3 == 0:
        pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
        old_pos = np.asarray(htf.old_positions(pos))
        traj_old.append(old_pos)
    integrator.step(step)
elapsed_time = (time.time() - initial_time) * unit.seconds
_logger.info(f'Done equilibrating, took: {elapsed_time / unit.seconds} seconds')

# Save traj
top_old = md.Topology.from_openmm(htf._topology_proposal.old_topology)
traj = md.Trajectory(np.array(traj_old), top_old)
traj.save(os.path.join(args.dir, f"{i}_{args.phase}_old.pdb"))
