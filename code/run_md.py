#!/bin/env python

import time
from simtk import openmm, unit
from simtk.openmm import app
from openmmtools.integrators import LangevinIntegrator
import argparse 
import mdtraj as md
from openmmforcefields.generators import SystemGenerator
from tqdm import tqdm
import os
import pickle
import numpy as np

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='complex or apo')
args = parser.parse_args()

#
# Global simulation parameters
#

pressure = 1.0 * unit.atmospheres
temperature = 310 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 2.0 * unit.femtoseconds
splitting = 'V R O R V'
nsteps = 1000 # 2 ps
niterations = 2500 # 5 ns

# Read in htf
dir_num = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{dir_num}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)

equilibrated_pdb_filename = 'equilibrated.pdb'

# Make integrator
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

# Minimize
print('Minimizing energy...')
context = openmm.Context(htf._topology_proposal.old_system, integrator)
context.setPositions(htf.old_positions(htf.hybrid_positions))
print('  initial : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))
openmm.LocalEnergyMinimizer.minimize(context)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))

# Equilibrate
print('Equilibrating...')
initial_time = time.time()
positions = []
for iteration in tqdm(range(niterations)):
    integrator.step(nsteps)
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    positions.append(pos.value_in_unit_system(unit.md_unit_system))

elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = niterations * nsteps * timestep
print('    Equilibration took %.3f s for %.3f ns (%8.3f ns/day)' % (elapsed_time / unit.seconds, simulation_time / unit.nanoseconds, simulation_time / elapsed_time * unit.day / unit.nanoseconds))

with open(os.path.join(args.dir, equilibrated_pdb_filename), 'w') as outfile:
    app.PDBFile.writeFile(htf._topology_proposal.old_topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))

# Save trajs
with open(os.path.join(args.dir, "positions.npy"), 'wb') as f:
    np.save(f, positions)
traj = md.Trajectory(np.array(np.array(positions)), md.Topology.from_openmm(htf._topology_proposal.old_topology))
traj.save(os.path.join(args.dir, "traj.dcd"))

