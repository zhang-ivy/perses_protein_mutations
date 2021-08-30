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
import pickle

# Read args
parser = argparse.ArgumentParser(description='analyze repex cache')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA, corresponding to the endstate')
parser.add_argument('endstate', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('length', type=int, help='5 or 10 ns')
args = parser.parse_args()

def get_dihedrals_all_replicas(i, aa, length, out_dir, topology, dihedral_indices):
    
    # From Hannah: https://github.com/hannahbrucemacdonald/endstate_pdbs/blob/master/scripts/input_for_pol_calc.py
    from perses.analysis.utils import open_netcdf
    nc = open_netcdf(os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}.nc"))
    nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_checkpoint.nc"))
    checkpoint_interval = nc_checkpoint.CheckpointInterval
    all_positions = nc_checkpoint.variables['positions']
    n_iter, n_replicas, n_atoms, _ = np.shape(all_positions)
    
    from tqdm import tqdm
    dihedrals_master = []
    for i in [0, 6, 11]:
        index = i # of replica
        all_pos = np.zeros(shape=(n_iter, topology.getNumAtoms(), 3))
        for iteration in tqdm(range(n_iter)):
            replica_id = np.where(nc.variables['states'][iteration*checkpoint_interval] == index)[0]
            pos = all_positions[iteration,replica_id,:,:][0] *unit.nanometers
            all_pos[iteration] = pos

        traj = md.Trajectory(np.array(all_pos), topology)
    #   dihedrals = np.sin(md.compute_dihedrals(traj, np.array([indices]))) 
        if len(dihedral_indices) == 4:
            print("computing thr dihedrals")
            dihedrals = md.compute_dihedrals(traj, np.array([dihedral_indices]))
        elif len(dihedral_indices) == 2:
            dihedrals_phi = md.compute_dihedrals(traj, np.array([dihedral_indices[0]]))
            dihedrals_psi = md.compute_dihedrals(traj, np.array([dihedral_indices[1]]))
            dihedrals = [dihedrals_phi, dihedrals_psi]
        dihedrals_master.append(dihedrals)
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

if args.name == 'THR':
    indices = [6, 7, 10, 12]
elif args.name == 'ALA':
    indices = [[6, 7, 8, 10], [14, 13, 8, 10]]

# Get topology
with open(os.path.join(args.dir, f"{i}_{args.name.lower()}_vanilla_topology.pickle"), "rb") as f:
    topology = pickle.load(f)    

dihedrals, n_iter = get_dihedrals_all_replicas(i, aa, length, out_dir, topology, indices)

for j, replica in tqdm(enumerate(dihedrals)):
    if len(replica) == (args.length*1000) + 1:
        plot_dihedrals(replica, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_correlated.png"))
        uncorrelated_old = plot_time_series(replica, n_iter, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_timeseries.png"))
        plot_dihedrals_uncorrelated(replica, uncorrelated_old, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}_decorrelated.png"))
    elif len(replica) ==2:
        plot_dihedrals(replica[0], os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}phi_correlated.png"))
        uncorrelated_old = plot_time_series(replica[0], n_iter, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}phi_timeseries.png"))
        plot_dihedrals_uncorrelated(replica[0], uncorrelated_old, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}phi_decorrelated.png"))
        plot_dihedrals(replica[1], os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}psi_correlated.png"))
        uncorrelated_old = plot_time_series(replica[1], n_iter, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}psi_timeseries.png"))
        plot_dihedrals_uncorrelated(replica[1], uncorrelated_old, os.path.join(out_dir, f"{i}_{args.phase}_{aa}_{length}_{j}_{aa}psi_decorrelated.png"))
