import argparse
import pickle
import os
from perses.annihilation.lambda_protocol import LambdaProtocol
import simtk.unit as unit
from openmmtools.multistate import MultiStateReporter
from perses.samplers.multistate import HybridRepexSampler
from openmmtools import mcmc

# Read WT and proposed residue abbreviations
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('input_file', type=str, help='path to hybrid topology factory')
args = parser.parse_args()

apo_htf = pickle.load(open(args.input_file, "rb" ))

# Build the hybrid repex samplers
suffix = 'run'; selection = 'not water'; checkpoint_interval = 10; n_states = 11; n_cycles = 5000
lambda_protocol = LambdaProtocol(functions='default')
reporter_file = f"../data/{os.path.basename(args.input_file)[:-3]}.nc"
reporter = MultiStateReporter(reporter_file, analysis_particle_indices = apo_htf.hybrid_topology.select(selection), checkpoint_interval = checkpoint_interval)
hss = HybridRepexSampler(mcmc_moves=mcmc.LangevinSplittingDynamicsMove(timestep= 4.0 * unit.femtoseconds,
                                                                      collision_rate=5.0 / unit.picosecond,
                                                                      n_steps=250,
                                                                      reassign_velocities=False,
                                                                      n_restart_attempts=20,
                                                                      splitting="V R R R O R R R V",
                                                                      constraint_tolerance=1e-06),
                                                                      hybrid_factory=apo_htf, online_analysis_interval=10)
hss.setup(n_states=n_states, temperature=300*unit.kelvin, storage_file=reporter, lambda_protocol=lambda_protocol, endstates=False)
hss.extend(n_cycles)

