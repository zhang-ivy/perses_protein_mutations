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
parser.add_argument('state', type=int, help="lambda value, i.e. 0 or 1")
args = parser.parse_args()

#
# Global simulation parameters
#

pressure = 1.0 * unit.atmospheres
temperature = 300 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 2.0 * unit.femtoseconds
# splitting="V R H O R V"
nsteps = 2500 # 5 ps
niterations = 4000 # 20 ns

# Read in htf
dir_num = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{dir_num}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
if args.state == 0:
	system = htf._topology_proposal.old_system
	positions = htf.old_positions(htf.hybrid_positions)
	topology = htf._topology_proposal.old_topology
elif args.state == 1:
	system = htf._topology_proposal.new_system
	positions = htf.new_positions(htf.hybrid_positions)
	topology = htf._topology_proposal.new_topology

equilibrated_pdb_filename = f'{args.state}_{args.phase}_equilibrated.pdb'

# Make integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Minimize
print('Minimizing energy...')
context = openmm.Context(system, integrator)
context.setPositions(positions)
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
    app.PDBFile.writeFile(topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))

# Save trajs
with open(os.path.join(args.dir, f"{args.state}_{args.phase}_positions.npy"), 'wb') as f:
    np.save(f, positions)
traj = md.Trajectory(np.array(np.array(positions)), md.Topology.from_openmm(topology))
traj.save(os.path.join(args.dir, f"{args.phase}_traj.dcd"))

