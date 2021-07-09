import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
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
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('length', type=float, help='neq switching time in ns')
parser.add_argument('is_rxn_field', type=bool, help='whether to use rxn field protocol')
args = parser.parse_args()

# Define lambda functions
_logger.info(f"Defining lambda fn")
x = 'lambda'
if not args.is_rxn_field:
    _logger.info(f"Using default PME protocol")

    ALCHEMICAL_FUNCTIONS = {
                             'lambda_sterics_core': x,
                             'lambda_electrostatics_core': x,
                             'lambda_sterics_insert': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_sterics_delete': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_insert': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_delete': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_bonds': x,
                             'lambda_angles': x,
                             'lambda_torsions': x}

    # inflection1, inflection2, inflection3, inflection4 = 0.2, 0.4, 0.6, 0.8
    # ALCHEMICAL_FUNCTIONS = {
    #                              'lambda_sterics_core': x,
    #                              'lambda_electrostatics_core': x,
    #                              'lambda_sterics_insert': f"select(step({x} - {inflection3}), select(step({x}-{inflection4}), 1, (1/{inflection1})*({x}-{inflection3})), 0.0)",
    #                              'lambda_sterics_delete': f"select(step({x} - {inflection1}), select(step({x} - {inflection2}), 1, (1/{inflection1})*({x}-{inflection1})), 0.0)",
    #                              'lambda_electrostatics_insert': f"select(step({x} - {inflection4}), (1/{inflection1}) * ({x} - {inflection4}), 0.0)",
    #                              'lambda_electrostatics_delete': f"select(step({x} - {inflection1}), 1.0, (1/{inflection1})*{x})",
    #                              'lambda_bonds': x,
    #                              'lambda_angles': x,
    #                              'lambda_torsions': x}
else:
    _logger.info(f"Using rxn field protocol")
    ALCHEMICAL_FUNCTIONS = {'lambda_0_bonds_old': f'abs(1 - {x})',
                             'lambda_0_bonds_new': f'abs({x})',
                             'lambda_0_angles_old': f'abs(1 - {x})',
                             'lambda_0_angles_new': f'abs({x})',
                             'lambda_0_torsions_old': f'abs(1 - {x})',
                             'lambda_0_torsions_new': f'abs({x})',
                             'lambda_0_electrostatics_old':f'abs(1 - {x})',
                             'lambda_0_electrostatics_new': f'abs({x})',
                             'lambda_0_electrostatics_exceptions_old': f'abs(1 - {x})',
                             'lambda_0_electrostatics_exceptions_new': f'abs({x})',
                             'lambda_0_sterics_old': f'abs(1 - {x})',
                             'lambda_0_sterics_new': f'abs({x})',
                             'lambda_0_sterics_exceptions_old': f'abs(1 - {x})',
                             'lambda_0_sterics_exceptions_new': f'abs({x})'}

# Define simulation parameters
nsteps_eq = 25000 # 100 ps
nsteps_neq = int(args.length*250000) # 1 ns
neq_splitting='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
temperature = 300.0 * unit.kelvin

# Read in vanilla htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
    positions = htf.hybrid_positions
    system = htf.hybrid_system
    box_vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()

# Set all heavy atom masses to be 0
heavy_atoms = []
for atom in htf.hybrid_topology.atoms:
    if atom.element.name != 'hydrogen' and atom.residue.name not in ['HOH', 'Na+', 'Cl-']:
        system.setParticleMass(atom.index, 0.0)
        heavy_atoms.append(atom.index)

for i in range(system.getNumConstraints() - 1, -1, -1):
    p1, p2, distance = system.getConstraintParameters(i)
    if p1 in heavy_atoms or p2 in heavy_atoms:
        system.removeConstraint(i)

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
forward_works = [integrator.get_protocol_work(dimensionless=True)]
for fwd_step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    _logger.info(f"Forward neq: {fwd_step*2500} completed")

    forward_works.append(integrator.get_protocol_work(dimensionless=True))
    
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    old_traj = md.Trajectory(old_pos, md.Topology.from_openmm(htf._topology_proposal.old_topology))
    old_pos_solute = old_traj.atom_slice(old_traj.top.select("not water")).xyz[0]
    
    new_pos = np.asarray(htf.new_positions(pos))
    new_traj = md.Trajectory(new_pos, md.Topology.from_openmm(htf._topology_proposal.new_topology))
    new_pos_solute = new_traj.atom_slice(new_traj.top.select("not water")).xyz[0]
    
    forward_neq_old.append(old_pos_solute)
    forward_neq_new.append(new_pos_solute)
forward_works_master.append(forward_works)

# Run eq reverse (1 -> 0)
integrator.step(nsteps_eq)

# Run neq reverse (1 -> 0)
reverse_works_master = list()
reverse_neq_old, reverse_neq_new = list(), list()
reverse_works = [integrator.get_protocol_work(dimensionless=True)]
for rev_step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    _logger.info(f"Reverse neq: {rev_step*2500} completed")
    
    reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    old_traj = md.Trajectory(old_pos, md.Topology.from_openmm(htf._topology_proposal.old_topology))
    old_pos_solute = old_traj.atom_slice(old_traj.top.select("not water")).xyz[0]
    
    new_pos = np.asarray(htf.new_positions(pos))
    new_traj = md.Trajectory(new_pos, md.Topology.from_openmm(htf._topology_proposal.new_topology))
    new_pos_solute = new_traj.atom_slice(new_traj.top.select("not water")).xyz[0]

    reverse_neq_old.append(old_pos_solute)
    reverse_neq_new.append(new_pos_solute)
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

