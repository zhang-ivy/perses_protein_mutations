import argparse
import pickle
import os
from perses.annihilation.lambda_protocol import LambdaProtocol
import simtk.unit as unit
from openmmtools.multistate import MultiStateReporter
from perses.samplers.multistate import HybridRepexSampler
from openmmtools import mcmc
import logging
import numpy as np

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='complex or apo')
args = parser.parse_args()

i = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), "rb" ))

# Build the hybrid repex samplers
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
selection = 'not water'; checkpoint_interval = 10; n_states = 32; n_cycles = 5000 
# lambda_protocol = LambdaProtocol(functions='default')

# Define multiphase protocol
inflection1, inflection2, inflection3, inflection4 = 0.2, 0.4, 0.6, 0.8
multiphase = {'lambda_sterics_core':
                         lambda x: x,
                         'lambda_electrostatics_core':
                         lambda x: x,
                         'lambda_sterics_insert':
                         lambda x: 0.0 if x < inflection3 else ((1.0/inflection1)*(x-inflection3) if x < inflection4 else 1.0),
                         'lambda_sterics_delete':
                         lambda x: 0.0 if x < inflection1 else ((1.0/inflection1)*(x-inflection1) if x < inflection2 else 1.0),
                         'lambda_electrostatics_insert':
                         lambda x: 0.0 if x < inflection4 else (1.0/inflection1)*(round(x-inflection4, 2)),
                         'lambda_electrostatics_delete':
                         lambda x: (1.0/inflection1)*x if x < inflection1 else 1.0,
                         'lambda_bonds':
                         lambda x: x,
                         'lambda_angles':
                         lambda x: x,
                         'lambda_torsions':
                         lambda x: x
                         }
lambda_protocol = LambdaProtocol(functions=multiphase)
lambda_schedule = list(np.linspace(0.0,0.4, 14)) +  [0.5] + list(np.linspace(0.6,0.7, 4)) + list(np.linspace(0.71, 0.8, 6)) + list(np.linspace(0.81, 1.0, 7))

reporter_file = os.path.join(args.dir, f"{i}_{args.phase}.nc")
reporter = MultiStateReporter(reporter_file, analysis_particle_indices = htf.hybrid_topology.select(selection), checkpoint_interval = checkpoint_interval)
hss = HybridRepexSampler(mcmc_moves=mcmc.LangevinSplittingDynamicsMove(timestep= 4.0 * unit.femtoseconds,
                                                                      collision_rate=5.0 / unit.picosecond,
                                                                      n_steps=250,
                                                                      reassign_velocities=False,
                                                                      n_restart_attempts=20,
                                                                      splitting="V R R R O R R R V",
                                                                      constraint_tolerance=1e-06),
                                                                      replica_mixing_scheme='swap-neighbors',
                                                                      hybrid_factory=htf, online_analysis_interval=10)
hss.setup(n_states=n_states, temperature=298*unit.kelvin, storage_file=reporter, lambda_schedule=lambda_schedule, lambda_protocol=lambda_protocol)
hss.extend(n_cycles)

