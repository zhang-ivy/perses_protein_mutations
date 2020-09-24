import argparse
import math
from simtk import unit
from openmmtools import testsystems, states, mcmc, multistate
import os
import tempfile
import pickle
from perses.annihilation.lambda_protocol import RelativeAlchemicalState
import mdtraj as md
import numpy as np
from simtk.unit.quantity import Quantity
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator


# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
args = parser.parse_args()

i = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), "rb" ))

# Generate system
pdb = app.PDBFile(os.path.join(args.dir, f"{args.name.lower()}_vacuum.pdb"))

forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
barostat = None
system_generator = SystemGenerator(forcefields=forcefield_files,
                               barostat=barostat,
                               forcefield_kwargs={'removeCMMotion': False,
                                                    'ewaldErrorTolerance': 1e-4,
                                                    'constraints' : app.HBonds,
                                                    'hydrogenMass' : 4 * unit.amus},
                                periodic_forcefield_kwargs=None,
                                small_molecule_forcefield='gaff-2.11',
                                nonperiodic_forcefield_kwargs = {'nonbondedMethod': app.NoCutoff}, 
                                   molecules=None, 
                                   cache=None)
system = system_generator.create_system(pdb.topology) 
positions = unit.quantity.Quantity(value = np.array([list(atom_pos) for atom_pos in pdb.positions.value_in_unit_system(unit.md_unit_system)]), unit = unit.nanometers)


# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

thermodynamic_states = [states.ThermodynamicState(system=system, temperature=T)
                        for T in temperatures]

# Set up sampler
move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=500)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_states,
                  sampler_states=states.SamplerState(positions),
                  storage=reporter)
simulation.run()