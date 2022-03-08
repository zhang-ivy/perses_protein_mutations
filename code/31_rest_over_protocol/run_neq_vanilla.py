import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from openmmtools.constants import kB
from simtk import unit, openmm
import argparse
import os
import time
import mdtraj as md
from tqdm import tqdm
from openeye import oechem

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('sim_number', type=int, help='number in job name - 1')
parser.add_argument('eq_length', type=float, help='eq switching time in ns')
parser.add_argument('neq_length', type=float, help='neq switching time in ns')
parser.add_argument('chain_A', type=int, help='(first) chain index for which to to add a virtual bond for complex phase')
parser.add_argument('chain_B', type=int, help='(second) chain index for which to to add a virtual bond for complex phase')
args = parser.parse_args()

# Define simulation parameters
nsteps_eq = int(args.eq_length*250000) # 1 ns
nsteps_neq = int(args.neq_length*250000) # 1 ns
neq_splitting='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
temperature = 300.0 * unit.kelvin

# Define lambda functions
x = 'lambda'

ALCHEMICAL_FUNCTIONS = {
'lambda_rest_bonds': "1",
'lambda_rest_angles': "1",
'lambda_rest_torsions':"1",
'lambda_rest_electrostatics': "1",
'lambda_rest_electrostatics_exceptions': "1",
'lambda_rest_sterics':"1",
'lambda_rest_sterics_exceptions': "1", 
                         'lambda_alchemical_bonds_old': f'1 - {x}',
                         'lambda_alchemical_bonds_new': x,
                         'lambda_alchemical_angles_old': f'1 - {x}',
                         'lambda_alchemical_angles_new': x,
                         'lambda_alchemical_torsions_old': f'1 - {x}',
                         'lambda_alchemical_torsions_new': x,
                         'lambda_alchemical_electrostatics_old': f'1 - {x}',
                         'lambda_alchemical_electrostatics_new': x,
                         'lambda_alchemical_electrostatics_exceptions_old': f'1 - {x}',
                         'lambda_alchemical_electrostatics_exceptions_new': x,
                         'lambda_alchemical_electrostatics_reciprocal': x,
                         'lambda_alchemical_sterics_old': f'1 - {x}',
                         'lambda_alchemical_sterics_new': x,
                         'lambda_alchemical_sterics_exceptions_old': f'1 - {x}',
                         'lambda_alchemical_sterics_exceptions_new': x
                         }

# Read in vanilla htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
positions = htf.hybrid_positions
box_vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()

# Make sure LRC is set correctly
system = htf.hybrid_system
system.getForce(5).setUseLongRangeCorrection(False)
system.getForce(8).setUseDispersionCorrection(True)
_logger.info(f"CustomNonbondedForce_sterics use LRC? {system.getForce(5).getUseLongRangeCorrection()}")
_logger.info(f"NonbondedForce_sterics use LRC? {system.getForce(8).getUseDispersionCorrection()}")

# Add virtual bond for complex phase
if args.phase == 'complex':
    chains = list(htf.hybrid_topology.chains)
    atom_A = list(chains[args.chain_A].atoms)[0]
    atom_B = list(chains[args.chain_B].atoms)[0]
    force = openmm.CustomBondForce('0')
    force.addBond(atom_A.index, atom_B.index, [])
    system.addForce(force)
    _logger.info(f"Added virtual bond between {atom_A} and {atom_B}")

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep, temperature=temperature)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*box_vectors)
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run eq forward (0 -> 1)
integrator.step(nsteps_eq)

# Run neq forward (0 -> 1)
forward_works_master = list()
forward_neq_old, forward_neq_new = list(), list()
forward_neq_old_waters, forward_neq_new_waters = list(), list()
forward_neq_box_vectors = list()
forward_works = [integrator.get_protocol_work(dimensionless=True)]
for fwd_step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    _logger.info(f"Forward neq: {fwd_step*2500} completed")

    forward_works.append(integrator.get_protocol_work(dimensionless=True))
    
    state = context.getState(getPositions=True, enforcePeriodicBox=False)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    pos = state.getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    old_traj = md.Trajectory(old_pos, md.Topology.from_openmm(htf._topology_proposal.old_topology))
    old_pos_solute = old_traj.atom_slice(old_traj.top.select("not water")).xyz[0]
    
    new_pos = np.asarray(htf.new_positions(pos))
    new_traj = md.Trajectory(new_pos, md.Topology.from_openmm(htf._topology_proposal.new_topology))
    new_pos_solute = new_traj.atom_slice(new_traj.top.select("not water")).xyz[0]
    
    forward_neq_box_vectors.append(box_vectors)
    forward_neq_old.append(old_pos_solute)
    forward_neq_new.append(new_pos_solute)
    if fwd_step == int(nsteps_neq / 2500) - 1:
        forward_neq_old_waters.append(old_pos)
        forward_neq_new_waters.append(new_pos)
forward_works_master.append(forward_works)

for k, v in context.getParameters().items():
    print(k, v)

# Run eq reverse (1 -> 0)
integrator.step(nsteps_eq)

# Run neq reverse (1 -> 0)
reverse_works_master = list()
reverse_neq_old, reverse_neq_new = list(), list()
reverse_neq_old_waters, reverse_neq_new_waters = list(), list()
reverse_neq_box_vectors = list()
reverse_works = [integrator.get_protocol_work(dimensionless=True)]
for rev_step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    _logger.info(f"Reverse neq: {rev_step*2500} completed")
    
    reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    
    state = context.getState(getPositions=True, enforcePeriodicBox=False)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    pos = state.getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    old_traj = md.Trajectory(old_pos, md.Topology.from_openmm(htf._topology_proposal.old_topology))
    old_pos_solute = old_traj.atom_slice(old_traj.top.select("not water")).xyz[0]
    
    new_pos = np.asarray(htf.new_positions(pos))
    new_traj = md.Trajectory(new_pos, md.Topology.from_openmm(htf._topology_proposal.new_topology))
    new_pos_solute = new_traj.atom_slice(new_traj.top.select("not water")).xyz[0]

    reverse_neq_box_vectors.append(box_vectors)
    reverse_neq_old.append(old_pos_solute)
    reverse_neq_new.append(new_pos_solute)
    if rev_step == int(nsteps_neq / 2500) - 1:
        reverse_neq_old_waters.append(old_pos)
        reverse_neq_new_waters.append(new_pos)
reverse_works_master.append(reverse_works)

# Save works
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward.npy"), 'wb') as f:
    np.save(f, forward_works_master)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse.npy"), 'wb') as f:
    np.save(f, reverse_works_master)

with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_old.npy"), 'wb') as f:
    np.save(f, forward_neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_new.npy"), 'wb') as f:
    np.save(f, forward_neq_new)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_old.npy"), 'wb') as f:
    np.save(f, reverse_neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_new.npy"), 'wb') as f:
    np.save(f, reverse_neq_new)

with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_old_waters.npy"), 'wb') as f:
    np.save(f, forward_neq_old_waters)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_new_waters.npy"), 'wb') as f:
    np.save(f, forward_neq_new_waters)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_old_waters.npy"), 'wb') as f:
    np.save(f, reverse_neq_old_waters)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_new_waters.npy"), 'wb') as f:
    np.save(f, reverse_neq_new_waters)

with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_forward_neq_box_vectors.npy"), 'wb') as f:
    np.save(f, forward_neq_box_vectors)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_reverse_neq_box_vectors.npy"), 'wb') as f:
    np.save(f, reverse_neq_box_vectors)
