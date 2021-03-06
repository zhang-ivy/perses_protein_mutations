import argparse
from openmmtools import multistate
import os

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('length', type=int, help='in ns')
args = parser.parse_args()

# Set up sampler
simulation = multistate.ReplicaExchangeSampler.from_storage(os.path.join(args.dir, "18_vacuum_thr_5ns.nc"))

# Run t-repex
simulation.extend(n_iterations=args.length*1000)