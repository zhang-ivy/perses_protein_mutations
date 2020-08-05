from openmmtools.integrators import AlchemicalNonequilibriumLangevinIntegrator
from simtk import unit
from simtk import openmm
from perses.annihilation.lambda_protocol import RelativeAlchemicalState
import mdtraj as md
import numpy as np
import logging
import pickle
import argparse
import os
import time
from perses.analysis.utils import open_netcdf
from perses.dispersed import feptasks


# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('sim_number', type=str, help='number in job name - 1')
args = parser.parse_args()

# Define lambda functions
x = 'lambda'
DEFAULT_ALCHEMICAL_FUNCTIONS = {
                             'lambda_sterics_core': f"1 - {x}",
                             'lambda_electrostatics_core': f"1 - {x}",
                             'lambda_sterics_insert': f"1 - select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_sterics_delete': f"1 - select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_insert': f"1 - select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_delete': f"1 - select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_bonds': f"1 - {x}",
                             'lambda_angles': f"1 - {x}",
                             'lambda_torsions': f"1 - {x}"}

# Define simulation parameters
nsteps_neq= 20000 # 80 ps
neq_splitting='V R H O R V'
timestep=4.0 * unit.femtosecond
temperature = 300*unit.kelvin
collision_rate = 90/unit.picosecond
platform_name = 'CUDA'

# Read in htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
new_top = md.Topology.from_openmm(htf._topology_proposal.new_topology)
system = htf.hybrid_system

# Read in cache
with open(os.path.join(args.dir, f"ser_pos_hybrid.npy"), 'rb') as f:
    ser_pos_hybrid = np.load(f)

# Read in indices of uncorrelated snapshots
with open(os.path.join(args.dir, f"ser_indices.npy"), 'rb') as f:
    ser_indices = np.load(f)

# Get equilbrium snapshot
positions = ser_pos_hybrid[ser_indices[int(args.sim_number)]]

# Set up integrator
integrator = AlchemicalNonequilibriumLangevinIntegrator(alchemical_functions=DEFAULT_ALCHEMICAL_FUNCTIONS, nsteps_neq=nsteps_neq, splitting=neq_splitting, timestep=timestep)

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
reverse_works_master = list()
reverse_neq_old, reverse_neq_new = list(), list()

# reverse (1 -> 0)
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
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse.npy"), 'wb') as f:
    np.save(f, reverse_works_master)

# Save trajs
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_old.npy"), 'wb') as f:
    np.save(f, reverse_neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_new.npy"), 'wb') as f:
    np.save(f, reverse_neq_new)
