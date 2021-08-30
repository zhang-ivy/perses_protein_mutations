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

# Read args
parser = argparse.ArgumentParser(description='analyze repex cache')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA, corresponding to the endstate')
parser.add_argument('endstate', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('length', type=int, help='5 or 10 ns')
args = parser.parse_args()

def get_dihedrals_all_replicas(i, aa, length, out_dir, htf, dihedral_indices_new, dihedral_indices_old):
    new_top = md.Topology.from_openmm(htf._topology_proposal.new_topology)
    old_top = md.Topology.from_openmm(htf._topology_proposal.old_topology)
    
    # From Hannah: https://github.com/hannahbrucemacdonald/endstate_pdbs/blob/master/scripts/input_for_pol_calc.py
    from perses.analysis.utils import open_netcdf
    nc = open_netcdf(os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}.nc"))
    nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_checkpoint.nc"))
    checkpoint_interval = nc_checkpoint.CheckpointInterval
    all_positions = nc_checkpoint.variables['positions']
    n_iter, n_replicas, n_atoms, _ = np.shape(all_positions)
    
    from tqdm import tqdm
    dihedrals_master = []
    for i in [0, 5, 10]:
        index = i # of replica
        all_pos_new = np.zeros(shape=(n_iter, new_top.n_atoms, 3))
        all_pos_old = np.zeros(shape=(n_iter, old_top.n_atoms, 3))
#         all_pos_hybrid = np.zeros(shape=(n_iter, n_atoms, 3))
        for iteration in tqdm(range(n_iter)):
            replica_id = np.where(nc.variables['states'][iteration*checkpoint_interval] == index)[0]
            pos = all_positions[iteration,replica_id,:,:][0] *unit.nanometers
            all_pos_new[iteration] = htf.new_positions(pos).value_in_unit_system(unit.md_unit_system) # Get new positions only
#             all_pos_hybrid[iteration] = pos.value_in_unit_system(unit.md_unit_system) # Get hybrid positions
            all_pos_old[iteration] = htf.old_positions(pos).value_in_unit_system(unit.md_unit_system)

        dihedrals_all = []
        # write for loop for this part -- old and new atoms
        for pos, top, indices in zip([all_pos_new, all_pos_old], [new_top, old_top], [dihedral_indices_new, dihedral_indices_old]):
            traj = md.Trajectory(np.array(pos), top)
    #         dihedrals = np.sin(md.compute_dihedrals(traj, np.array([indices]))) 
            dihedrals = md.compute_dihedrals(traj, np.array([indices]))
            dihedrals_all.append(dihedrals)
        dihedrals_master.append(dihedrals_all)
    return dihedrals_master, n_iter
    
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

i = os.path.basename(os.path.dirname(args.dir))
aa = args.name.lower()
length = f'{args.length}ns'
out_dir = args.dir

with open(os.path.join(out_dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)

if args.name == 'THR':
    indices_old = [12, 10, 7, 6]
    indices_new = [13, 10, 7, 6]
elif args.name == 'SER':
    indices_old = [11, 10, 7, 6]
    indices_new = [15, 10, 7, 6]
dihedrals, n_iter = get_dihedrals_all_replicas(i, aa, length, out_dir, htf, indices_new, indices_old)

for j, replica in tqdm(enumerate(dihedrals)):
    dihedrals_new = replica[0]
    dihedrals_old = replica[1]
    plot_dihedrals(dihedrals_old, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_correlated.png"))
    uncorrelated_old = plot_time_series(dihedrals_old, n_iter, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_timeseries.png"))
    plot_dihedrals_uncorrelated(dihedrals_old, uncorrelated_old, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_decorrelated.png"))
    plot_dihedrals(dihedrals_new, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_ala_correlated.png"))
    uncorrelated_new = plot_time_series(dihedrals_new, n_iter, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_ala_timeseries.png"))
    plot_dihedrals_uncorrelated(dihedrals_new, uncorrelated_new, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_ala_decorrelated.png"))
