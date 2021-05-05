import argparse
from openmmtools.multistate import MultiStateReporter
import logging
import os
import numpy as np

# Read args
parser = argparse.ArgumentParser(description='resume rest')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('resid', type=str, help='resid of residue to be mutated')
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('sim_number', type=int, help='index of job array, starts at 1')
parser.add_argument('--total_iterations', type=int, help='total number of iterations desired')
args = parser.parse_args()

if args.sim_number == 1:
    phase = 'apo'
    name = args.old_aa_name
elif args.sim_number == 2:
    phase = 'apo'
    name = args.new_aa_name
elif args.sim_number == 3:
    phase = 'complex'
    name = args.old_aa_name
elif args.sim_number == 4:
    phase = 'complex'
    name = args.new_aa_name

length = 1

from openmmtools.multistate import ReplicaExchangeSampler
import mpiplus
class ReplicaExchangeSampler2(ReplicaExchangeSampler):
    @mpiplus.on_single_node(rank=0, broadcast_result=False, sync_nodes=False)
    @mpiplus.delayed_termination
    def _report_iteration_items(self):
        """
        Sub-function of :func:`_report_iteration` which handles all the actual individual item reporting in a
        sub-class friendly way. The final actions of writing timestamp, last-good-iteration, and syncing
        should be left to the :func:`_report_iteration` and subclasses should extend this function instead
        """
        replica_id = np.where(self._replica_thermodynamic_states == 0)[0][0]
        print("ITERATION: ", self._iteration)
        print("REPLICA THERMOSTATES ", self._replica_thermodynamic_states, type(self._replica_thermodynamic_states))
        print("REPLICA ID ", replica_id, type(replica_id))
        print("REPLICA SAMPLER STATES ", self._sampler_states)
        self._reporter.write_sampler_states([self._sampler_states[replica_id]], self._iteration)
        
        self._reporter.write_replica_thermodynamic_states(self._replica_thermodynamic_states, self._iteration)
        self._reporter.write_mcmc_moves(self._mcmc_moves)  # MCMCMoves can store internal statistics.
        self._reporter.write_energies(self._energy_thermodynamic_states, self._neighborhoods, self._energy_unsampled_states,
                                      self._iteration)
        self._reporter.write_mixing_statistics(self._n_accepted_matrix, self._n_proposed_matrix, self._iteration)

# Load repex simulation
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)
i = os.path.basename(os.path.dirname(args.dir))
reporter_file = os.path.join(args.dir, f"{i}_{phase}_{name.lower()}_{length}ns.nc")
reporter = MultiStateReporter(reporter_file, checkpoint_interval=10)
simulation = ReplicaExchangeSampler2.from_storage(reporter)

# Determine how many more iterations are needed
total_iterations = args.total_iterations if args.total_iterations else 1000
iterations =  total_iterations - simulation.iteration

# Resume simulation
simulation.extend(n_iterations=iterations)
