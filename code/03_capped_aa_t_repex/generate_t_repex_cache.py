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

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

alchemical_state = RelativeAlchemicalState.from_system(htf.hybrid_system)
alchemical_state.set_alchemical_parameters(args.state)

thermodynamic_states = [states.ThermodynamicState(system=htf.hybrid_system, temperature=T)
                        for T in temperatures]
compound_states = [states.CompoundThermodynamicState(thermodynamic_state=thermodynamic_state, composable_states=[alchemical_state])
                        for thermodynamic_state in thermodynamic_states]

# Create sampler state
sampler_state = states.SamplerState(htf.hybrid_positions)
sampler_state.box_vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()

# Set up sampler
move = mcmc.GHMCMove(timestep=4.0*unit.femtoseconds, n_steps=250)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=sampler_state,
                  storage=reporter)
simulation.run()