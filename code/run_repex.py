import argparse
import pickle
import os
from perses.annihilation.lambda_protocol import LambdaProtocol
import simtk.unit as unit
from openmmtools.multistate import MultiStateReporter
from perses.samplers.multistate import HybridRepexSampler
from openmmtools import mcmc
import logging
import datetime

# Set up logger
class TimeFilter(logging.Filter):
    def filter(self, record):
        try:
            last = self.last
        except AttributeError:
            last = record.relativeCreated
        delta = datetime.datetime.fromtimestamp(record.relativeCreated/1000.0) - datetime.datetime.fromtimestamp(last/1000.0)
        record.relative = '{0:.2f}'.format(delta.seconds + delta.microseconds/1000000.0)
        self.last = record.relativeCreated
        return True
fmt = logging.Formatter(fmt="%(asctime)s:(%(relative)ss):%(name)s:%(message)s")
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
[hndl.addFilter(TimeFilter()) for hndl in _logger.handlers]
[hndl.setFormatter(fmt) for hndl in _logger.handlers]

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
args = parser.parse_args()

simulation_number = os.path.basename(os.path.dirname(args.dir))
apo_htf = pickle.load(open(os.path.join(args.dir, f"{simulation_number}_{args.phase}.pickle"), "rb" ))

# Build the hybrid repex samplers
suffix = 'run'; selection = 'not water'; checkpoint_interval = 10; n_states = 11; n_cycles = 5000
lambda_protocol = LambdaProtocol(functions='default')
reporter_file = os.path.join(args.dir, f"{simulation_number}_{args.phase}.nc")
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

