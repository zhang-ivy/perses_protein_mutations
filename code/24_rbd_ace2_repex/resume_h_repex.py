import argparse
from perses.samplers.multistate import HybridRepexSampler
from openmmtools.multistate import MultiStateReporter
import logging
import os

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='complex or apo')
parser.add_argument('--total_iterations', type=int, help='total number of iterations desired')
args = parser.parse_args()

# Load repex simulation
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
i = os.path.basename(os.path.dirname(args.dir))
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}.nc")
reporter = MultiStateReporter(reporter_file, checkpoint_interval=10)
simulation = HybridRepexSampler.from_storage(reporter)

# Determine how many more iterations are needed
total_iterations = args.total_iterations if args.total_iterations else 5000
iterations =  total_iterations - simulation.iteration

# Resume simulation
simulation.extend(n_iterations=iterations)
