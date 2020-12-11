from simtk import unit
import os
import pickle
import mdtraj as md
import numpy as np
import logging 
import argparse

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA, corresponding to the endstate')
parser.add_argument('endstate', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('index', type=int, help='index of the repex state to grab traj for')
args = parser.parse_args()

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

def get_trajs_for_state(i, aa, phase, length, out_dir, htf, index):
    new_top = md.Topology.from_openmm(htf._topology_proposal.new_topology)
    old_top = md.Topology.from_openmm(htf._topology_proposal.old_topology)
    
    # From Hannah: https://github.com/hannahbrucemacdonald/endstate_pdbs/blob/master/scripts/input_for_pol_calc.py
    from perses.analysis.utils import open_netcdf
    nc = open_netcdf(os.path.join(out_dir, f"{i}_{phase}_{aa}_{length}.nc"))
    nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{i}_{phase}_{aa}_{length}_checkpoint.nc"))
    checkpoint_interval = nc_checkpoint.CheckpointInterval
    all_positions = nc_checkpoint.variables['positions']
    n_iter, n_replicas, n_atoms, _ = np.shape(all_positions)
    
    from tqdm import tqdm
    all_pos_new = np.zeros(shape=(n_iter, new_top.n_atoms, 3))
    all_pos_old = np.zeros(shape=(n_iter, old_top.n_atoms, 3))
    for iteration in tqdm(range(n_iter)):
        replica_id = np.where(nc.variables['states'][iteration*checkpoint_interval] == index)[0]
        pos = all_positions[iteration,replica_id,:,:][0] *unit.nanometers
        all_pos_new[iteration] = htf.new_positions(pos).value_in_unit_system(unit.md_unit_system) # Get new positions only
        all_pos_old[iteration] = htf.old_positions(pos).value_in_unit_system(unit.md_unit_system)

    return all_pos_new, all_pos_old

i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.endstate}.pickle"), "rb") as f:
    htf = pickle.load(f)

new, old = get_trajs_for_state(i, args.name.lower(), args.phase, args.length, args.dir, htf, args.index)

traj_new = md.Trajectory(new, md.Topology.from_openmm(htf._topology_proposal.new_topology))
traj_new.save(os.path.join(args.dir, f"new_{args.index}.dcd"))
traj_new[0].save(os.path.join(args.dir, f"new_{args.index}.pdb"))

traj_old = md.Trajectory(old, md.Topology.from_openmm(htf._topology_proposal.old_topology))
traj_old.save(os.path.join(args.dir, f"old_{args.index}.dcd"))
traj_old[0].save(os.path.join(args.dir, f"old_{args.index}.pdb"))