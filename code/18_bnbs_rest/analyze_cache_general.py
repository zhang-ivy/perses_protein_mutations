import math
from simtk import unit
import os
import tempfile
import pickle
import mdtraj as md
import numpy as np
from simtk.unit.quantity import Quantity
import logging 

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

from matplotlib import pyplot as plt
from simtk.openmm import app
from tqdm import tqdm
import argparse
import random

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('resid', type=str, help='resid of residue to be mutated')
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('sim_number', type=int, help='index of job array, starts at 1')
args = parser.parse_args()

if args.sim_number == 1:
    phase = 'apo'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 2:
    phase = 'apo'
    name = args.new_aa_name
    state = 1
elif args.sim_number == 3:
    phase = 'complex'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 4:
    phase = 'complex'
    name = args.new_aa_name
    state = 1

length = 1
i = os.path.basename(os.path.dirname(args.dir))
out_dir = args.dir

def get_dihedrals(i, name, length, out_dir, htf, dihedral_indices_new, dihedral_indices_old):
    new_top = md.Topology.from_openmm(htf._topology_proposal.new_topology)
    old_top = md.Topology.from_openmm(htf._topology_proposal.old_topology)
    
    # From Hannah: https://github.com/hannahbrucemacdonald/endstate_pdbs/blob/master/scripts/input_for_pol_calc.py
    from perses.analysis.utils import open_netcdf
    nc = open_netcdf(os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns.nc"))
    nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_checkpoint.nc"))
    checkpoint_interval = nc_checkpoint.CheckpointInterval
    all_positions = nc_checkpoint.variables['positions']
    n_iter, n_replicas, n_atoms, _ = np.shape(all_positions)
    
    from tqdm import tqdm
    all_pos_new = np.zeros(shape=(n_iter, new_top.n_atoms, 3))
    all_pos_old = np.zeros(shape=(n_iter, old_top.n_atoms, 3))
    all_pos_hybrid = np.zeros(shape=(n_iter, n_atoms, 3))
    for iteration in tqdm(range(n_iter)):
        replica_id = np.where(nc.variables['states'][iteration*checkpoint_interval] == 0)[0]
        pos = all_positions[iteration,replica_id,:,:][0] *unit.nanometers
        all_pos_new[iteration] = htf.new_positions(pos).value_in_unit_system(unit.md_unit_system) # Get new positions only
        all_pos_hybrid[iteration] = pos # Get hybrid positions
        all_pos_old[iteration] = htf.old_positions(pos).value_in_unit_system(unit.md_unit_system)

    dihedrals_all = []
    for pos, top, indices in zip([all_pos_new, all_pos_old], [new_top, old_top], [dihedral_indices_new, dihedral_indices_old]):
        traj = md.Trajectory(np.array(pos), top)
        dihedrals = md.compute_dihedrals(traj, np.array([indices]))
        dihedrals_all.append(dihedrals)
    
    return dihedrals_all, n_iter, all_pos_hybrid
    
def plot_dihedrals(dihedrals, outfile):
    # Plot histogram with error bars : https://stackoverflow.com/questions/35390276/how-to-add-error-bars-to-histogram-diagram-in-python
    entries, edges, _ = plt.hist(dihedrals)
    bin_centers = 0.5 * (edges[:-1] + edges[1:]) # calculate bin centers
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='r.') # draw errobars, use the sqrt error
    plt.xlim(-np.pi, np.pi)
    plt.savefig(outfile, dpi=300)
    plt.close()
    
def plot_time_series(dihedrals, n_iter, outfile):
    from perses.dispersed import feptasks
    t0, g, neff_max, a_t, uncorrelated_indices = feptasks.compute_timeseries(dihedrals)

    plt.scatter(range(n_iter), dihedrals)
    plt.ylabel("dihedral")
    plt.xlabel("iteration number")
    plt.ylim(-np.pi, np.pi)
    plt.savefig(outfile, dpi=300)
    plt.close()
    
    return uncorrelated_indices
    
def plot_dihedrals_uncorrelated(dihedrals, uncorrelated_indices, outfile):
    # Plot histogram with error bars : https://stackoverflow.com/questions/35390276/how-to-add-error-bars-to-histogram-diagram-in-python
    entries, edges, _ = plt.hist(dihedrals[uncorrelated_indices])
    bin_centers = 0.5 * (edges[:-1] + edges[1:]) # calculate bin centers
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='r.') # draw errobars, use the sqrt error
    plt.xlim(-np.pi, np.pi)
    plt.savefig(outfile, dpi=300)
    plt.close()


with open(os.path.join(out_dir, f"{i}_{phase}_{state}.pickle"), 'rb') as f:
    htf = pickle.load(f)

thr_dihedral = ['N', 'CA', 'CB', 'OG']
ring_dihedral = ['N', 'CA', 'CB', 'CG']
ala_dihedral = ['N', 'CA', 'CB', 'HB1']

dihedral_atoms = []
for aa_name in [args.old_aa_name, args.new_aa_name]:
    if aa_name in ["PHE", "TYR", "TRP"]:
        dihedral_atoms.append(ring_dihedral)
    elif aa_name == "THR":
        dihedral_atoms.append(thr_dihedral)
    elif aa_name == "ALA":
        dihedral_atoms.append(ala_dihedral)

indices_old = [atom.index for atom in list(htf._topology_proposal.old_topology.residues())[int(args.resid)].atoms() if atom.name in dihedral_atoms[0]]
indices_new = [atom.index for atom in list(htf._topology_proposal.new_topology.residues())[int(args.resid)].atoms() if atom.name in dihedral_atoms[1]]
_logger.info(f"old indices: {indices_old}")
_logger.info(f"new indices: {indices_new}")

dihedrals, n_iter, all_pos_hybrid = get_dihedrals(i, name, length, out_dir, htf, indices_new, indices_old)

# Plot 
dihedrals_new = dihedrals[0]
dihedrals_old = dihedrals[1]
plot_dihedrals(dihedrals_old, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.old_aa_name.lower()}_correlated.png"))
uncorrelated_old = plot_time_series(dihedrals_old, n_iter, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.old_aa_name.lower()}_timeseries.png"))
plot_dihedrals_uncorrelated(dihedrals_old, uncorrelated_old, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.old_aa_name.lower()}_decorrelated.png"))
plot_dihedrals(dihedrals_new, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.new_aa_name.lower()}_correlated.png"))
uncorrelated_new = plot_time_series(dihedrals_new, n_iter, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.new_aa_name.lower()}_timeseries.png"))
plot_dihedrals_uncorrelated(dihedrals_new, uncorrelated_new, os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_{args.new_aa_name.lower()}_decorrelated.png"))

# Save 100 random uncorrelated hybrid pos snapshots
if name == args.new_aa_name:
    uncorrelated_indices = uncorrelated_new
elif name == args.old_aa_name:
    uncorrelated_indices = uncorrelated_old
else:
    raise Exception("Your specified amino acid did not match the old or new aa names")
subset_indices = random.choices(uncorrelated_indices, k=100) # Choose 100 random indices from uncorrelated indices
_logger.info(f"randomly chosen indices: {subset_indices}")
subset_pos = all_pos_hybrid[subset_indices] # Make array of hybrid positions for 100 uncorrelated indices
with open(os.path.join(out_dir, f"{i}_{phase}_{name.lower()}_{length}ns_snapshots.npy"), 'wb') as f:
    np.save(f, subset_pos)
