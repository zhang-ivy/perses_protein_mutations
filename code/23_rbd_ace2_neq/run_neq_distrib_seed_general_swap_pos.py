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
import copy

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
parser.add_argument('--cache', type=int, default=1, help='length of rest cache in ns')
args = parser.parse_args()

# Define lambda functions
x = 'lambda'
# DEFAULT_ALCHEMICAL_FUNCTIONS = {
#                              'lambda_sterics_core': x,
#                              'lambda_electrostatics_core': x,
#                              'lambda_sterics_insert': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
#                              'lambda_sterics_delete': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
#                              'lambda_electrostatics_insert': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
#                              'lambda_electrostatics_delete': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
#                              'lambda_bonds': x,
#                              'lambda_angles': x,
#                              'lambda_torsions': x}

inflection1, inflection2, inflection3, inflection4 = 0.2, 0.4, 0.6, 0.8
DEFAULT_ALCHEMICAL_FUNCTIONS = {
                             'lambda_sterics_core': x,
                             'lambda_electrostatics_core': x,
                             'lambda_sterics_insert': f"select(step({x} - {inflection3}), select(step({x}-{inflection4}), 1, (1/{inflection1})*({x}-{inflection3})), 0.0)",
                             'lambda_sterics_delete': f"select(step({x} - {inflection1}), select(step({x} - {inflection2}), 1, (1/{inflection1})*({x}-{inflection1})), 0.0)",
                             'lambda_electrostatics_insert': f"select(step({x} - {inflection4}), (1/{inflection1}) * ({x} - {inflection4}), 0.0)",
                             'lambda_electrostatics_delete': f"select(step({x} - {inflection1}), 1.0, (1/{inflection1})*{x})",
                             'lambda_bonds': x,
                             'lambda_angles': x,
                             'lambda_torsions': x}

# Define simulation parameters
nsteps_eq = 10
nsteps_neq = int(args.length*250000) # 1 ns
neq_splitting='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
cache_length = args.cache if args.cache else 1 
temperature = 298.0 * unit.kelvin

# Read in vanilla htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)

# Read in lambda = 0 cache
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.old_aa_name}_{cache_length}ns_snapshots.npy"), 'rb') as f:
    subset_pos = np.load(f)
positions = subset_pos[args.sim_number]
system = htf.hybrid_system

# Swap positions of water to be turned into ion with a random water
def choose_water_atoms(charge_diff, old_indices, hybrid_positions, hybrid_topology, radius=0.8):
    import mdtraj as md
    # Create trajectory
    traj = md.Trajectory(hybrid_positions, hybrid_topology)
    water_atoms = traj.topology.select(f"water")
    query_atoms = traj.top.select('protein')

    while True:
        # Get water atoms within radius of protein
        neighboring_atoms = md.compute_neighbors(traj, radius, query_atoms, haystack_indices=water_atoms)[0]

        # Get water atoms outside of radius of protein
        nonneighboring_residues = set([atom.residue.index for atom in traj.topology.atoms if (atom.index in water_atoms) and (atom.index not in neighboring_atoms)])
        assert len(nonneighboring_residues) > 0, "there are no available nonneighboring waters"
        # Choose N random nonneighboring waters, where N is determined based on the charge_diff
        choice_residues = np.random.choice(list(nonneighboring_residues), size=abs(charge_diff), replace=False)

        # Get the atom indices in the water(s)
        choice_indices = np.array([[atom.index for atom in traj.topology.residue(res).atoms] for res in choice_residues])
        
        if not np.array_equal(choice_indices, np.array(old_indices)):
            break
        
    return np.ndarray.flatten(choice_indices)
og_indices = [150813, 150814, 150815] if args.phase == 'complex' else [37525, 37526, 37527]
swap_indices = choose_water_atoms(-1, og_indices, htf.hybrid_positions, htf.hybrid_topology)
for og, swap in zip(og_indices, swap_indices):
    temp = copy.copy(positions[og])
    positions[og] = positions[swap]
    positions[swap] = temp

# Read in lambda = 0 cache box vectors
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.old_aa_name}_{cache_length}ns_box_vectors.npy"), 'rb') as f:
    subset_box_vectors = np.load(f)
box_vectors = subset_box_vectors[args.sim_number][0]

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(DEFAULT_ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep, temperature=temperature)

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

# Read in lambda = 1 cache, if necessary
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.new_aa_name}_{cache_length}ns_snapshots.npy"), 'rb') as f:
    subset_pos = np.load(f)
positions = subset_pos[args.sim_number]

# Swap positions of water to be turned into ion with a random water
for og, swap in zip(og_indices, swap_indices):
    temp = copy.copy(positions[og])
    positions[og] = positions[swap]
    positions[swap] = temp

# Read in lambda = 1 cache box vectors
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.new_aa_name}_{cache_length}ns_box_vectors.npy"), 'rb') as f:
    subset_box_vectors = np.load(f)
box_vectors = subset_box_vectors[args.sim_number][0]

context.setPeriodicBoxVectors(*box_vectors)
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

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

