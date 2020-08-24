import argparse
from openmmtools import testsystems, states, mcmc, multistate
import os

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
args = parser.parse_args()

# Set up sampler
simulation = ReplicaExchangeSampler.from_storage(os.path.join(args.dir, "18_vacuum_thr_5ns.nc"))

# Run t-repex
simulation.extend(n_iterations=5000)