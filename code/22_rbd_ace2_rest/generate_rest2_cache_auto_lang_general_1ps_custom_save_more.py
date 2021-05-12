import pickle
import os
from perses.annihilation.rest import RESTTopologyFactory
from perses.annihilation.lambda_protocol import RESTState
from openmmtools.states import SamplerState, ThermodynamicState, CompoundThermodynamicState
from openmmtools import cache, utils
from perses.dispersed.utils import configure_platform
cache.global_context_cache.platform = configure_platform(utils.get_fastest_platform().getName())
from simtk import openmm, unit
import math
from openmmtools.constants import kB
from openmmtools import mcmc, multistate
import argparse
import copy
from perses.dispersed import feptasks
import mdtraj as md
import numpy as np
from perses.app.relative_point_mutation_setup import PointMutationExecutor

# Set up logger
import logging
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('resid', type=str, help='resid of residue to be mutated')
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('sim_number', type=int, help='index of job array, starts at 1')
parser.add_argument('t_max', type=int, help='max temp for rest')
parser.add_argument('radius', type=float, help='radius in nm for rest region')
args = parser.parse_args()

if args.sim_number == 1:
    phase = 'apo'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 2:
    phase = 'apo'
    name = args.new_aa_name
    state = 1
elif args.sim_number == 3:
    phase = 'complex'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 4:
    phase = 'complex'
    name = args.new_aa_name
    state = 1

length = 1
move_length = 1
timestep = 4
radius = args.radius

# Load rhtf 
i = os.path.basename(os.path.dirname(args.dir))
path = os.path.join(args.dir, f"{i}_{phase}_{state}.pickle")
_logger.info(f"path: {path}")
htf = pickle.load(open(path, "rb" ))

# Build REST factory
_logger.info("Generating REST factory")
_logger.info(f"radius:{radius} nm")
# query_indices = [atom.index for atom in list(htf.hybrid_topology.residues)[int(args.resid)].atoms]
for res in htf.hybrid_topology.residues:
    if res.resSeq == int(args.resid) and res.chain.index == 0:
        mutated_res = res
query_indices = [atom.index for atom in mutated_res.atoms]
_logger.info(f"query indices {query_indices}")
traj = md.Trajectory(np.array(htf.hybrid_positions), htf.hybrid_topology)
solute_atoms = list(traj.topology.select("is_protein"))
rest_atoms = list(md.compute_neighbors(traj, radius, query_indices, haystack_indices=solute_atoms)[0])
_logger.info(f"rest atoms {rest_atoms}")
factory = RESTTopologyFactory(htf.hybrid_system, solute_region=rest_atoms)

_logger.info("Generating REST states")
# Get REST system
REST_system = factory.REST_system

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = args.t_max * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

# Create reference thermodynamic state
lambda_zero_alchemical_state = RESTState.from_system(REST_system)
thermostate = ThermodynamicState(REST_system, temperature=T_min)
compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[lambda_zero_alchemical_state])

context_cache = cache.ContextCache()

# Create thermodynamics states
sampler_state =  SamplerState(htf.hybrid_positions, box_vectors=htf.hybrid_system.getDefaultPeriodicBoxVectors())
beta_0 = 1/(kB*T_min)
thermodynamic_state_list = []
sampler_state_list = []
for temperature in temperatures:
    beta_m = 1/(kB*temperature)
    compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
    compound_thermodynamic_state_copy.set_alchemical_parameters(beta_0, beta_m)
    thermodynamic_state_list.append(compound_thermodynamic_state_copy)

    # now generating a sampler_state for each thermodynamic state, with relaxed positions
    # context, context_integrator = context_cache.get_context(compound_thermodynamic_state_copy)
    feptasks.minimize(compound_thermodynamic_state_copy, sampler_state, max_iterations=0)
    sampler_state_list.append(copy.deepcopy(sampler_state))

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
        replica_id_0 = np.where(self._replica_thermodynamic_states == 0)[0][0]
        replica_id_23 = np.where(self._replica_thermodynamic_states == 23)[0][0]
        print("ITERATION: ", self._iteration)
        print("REPLICA THERMOSTATES ", self._replica_thermodynamic_states, type(self._replica_thermodynamic_states))
        print("REPLICA ID ", replica_id, type(replica_id_23))
        self._reporter.write_sampler_states([self._sampler_states[replica_id_0], self._sampler_states[replica_id_23]], self._iteration)
        
        self._reporter.write_replica_thermodynamic_states(self._replica_thermodynamic_states, self._iteration)
        self._reporter.write_mcmc_moves(self._mcmc_moves)  # MCMCMoves can store internal statistics.
        self._reporter.write_energies(self._energy_thermodynamic_states, self._neighborhoods, self._energy_unsampled_states,
                                      self._iteration)
        self._reporter.write_mixing_statistics(self._n_accepted_matrix, self._n_proposed_matrix, self._iteration)

# Set up sampler
_logger.setLevel(logging.DEBUG)
_logger.info("About to start repex")
print(f"move steps: {int((move_length*1000)/timestep)}")
print(f"timestep: {timestep} fs")
move = mcmc.LangevinSplittingDynamicsMove(timestep=timestep*unit.femtoseconds, n_steps=int((move_length*1000)/timestep))
simulation = ReplicaExchangeSampler2(mcmc_moves=move, number_of_iterations=length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{phase}_{name.lower()}_{length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()
