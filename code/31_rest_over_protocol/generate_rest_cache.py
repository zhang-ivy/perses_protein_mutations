import pickle
import os
import argparse
import copy
import mdtraj as md
import numpy as np
import math

from perses.annihilation.lambda_protocol import RESTCapableRelativeAlchemicalState, RESTCapableLambdaProtocol
from perses.dispersed import feptasks
from perses.dispersed.utils import configure_platform

from openmmtools.states import SamplerState, ThermodynamicState, CompoundThermodynamicState
from openmmtools import cache, utils
from openmmtools.constants import kB
from openmmtools import mcmc, multistate
cache.global_context_cache.platform = configure_platform(utils.get_fastest_platform().getName())
context_cache = cache.ContextCache(capacity=None, time_to_live=None)

from simtk.openmm import unit, app
from simtk import openmm

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
parser.add_argument('n_replicas', type=int, help='number of replicas')
parser.add_argument('length', type=int, help='length in ns for each replica')
parser.add_argument('chain_A', type=int, help='(first) chain index for which to to add a virtual bond for complex phase')
parser.add_argument('chain_B', type=int, help='(second) chain index for which to to add a virtual bond for complex phase')
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

length = args.length
move_length = 1
timestep = 4

# Load rest htf 
i = os.path.basename(os.path.dirname(args.dir))
path = os.path.join(args.dir, f"{i}_{phase}.pickle")
_logger.info(f"Loading htf at path: {path}")
htf = pickle.load(open(path, "rb" ))
positions = htf.hybrid_positions
REST_system = htf.hybrid_system
_logger.info(f"REST region has {len(htf._rest_region)} atoms")

# Make sure LRC is set correctly
REST_system.getForce(5).setUseLongRangeCorrection(False)
REST_system.getForce(8).setUseDispersionCorrection(True)
_logger.info(f"CustomNonbondedForce_sterics use LRC? {REST_system.getForce(5).getUseLongRangeCorrection()}")
_logger.info(f"NonbondedForce_sterics use LRC? {REST_system.getForce(8).getUseDispersionCorrection()}")

# Add virtual bond for complex phase
if phase == 'complex':
    chains = list(htf.hybrid_topology.chains)
    atom_A = list(chains[args.chain_A].atoms)[0]
    atom_B = list(chains[args.chain_B].atoms)[0]
    force = openmm.CustomBondForce('0')
    force.addBond(atom_A.index, atom_B.index, [])
    REST_system.addForce(force)
    _logger.info(f"Added virtual bond between {atom_A} and {atom_B}")

# Generate temperatures for each state
_logger.info("Generating temperatures")
n_replicas = args.n_replicas  # Number of temperature replicas.
T_min = 300.0 * unit.kelvin  # Minimum temperature.
T_max = args.t_max * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

# Create reference thermodynamic state
_logger.info("Generating reference thermostate")
lambda_protocol = RESTCapableLambdaProtocol(functions='no-alchemy')
lambda_zero_alchemical_state = RESTCapableRelativeAlchemicalState.from_system(REST_system)
thermostate = ThermodynamicState(REST_system, temperature=T_min)
compound_thermodynamic_state = CompoundThermodynamicState(thermostate,
                                                          composable_states=[lambda_zero_alchemical_state])

# Set alchemical parameters
_logger.info("Setting reference state alchemical parameters")
beta_0 = 1 / (kB * T_min)
beta_m = 1 / (kB * T_max)
global_lambda = 1
compound_thermodynamic_state.set_alchemical_parameters(global_lambda, beta_0, beta_m, lambda_protocol=lambda_protocol, endstate=state)
    
# Create thermodynamics states
_logger.info("Generating thermostates")
sampler_state =  SamplerState(positions, box_vectors=REST_system.getDefaultPeriodicBoxVectors())
thermodynamic_state_list = []
sampler_state_list = []
for temperature in temperatures:
    beta_m = 1/(kB*temperature)
    compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
    compound_thermodynamic_state_copy.set_alchemical_parameters(global_lambda, beta_0, beta_m, lambda_protocol=lambda_protocol, endstate=state)
    thermodynamic_state_list.append(compound_thermodynamic_state_copy)

    # now generating a sampler_state for each thermodynamic state, with relaxed positions
    _logger.info(f"Minimizing state for temp: {temperature}")
    _logger.info(f"State alchemical parameters: {compound_thermodynamic_state_copy.lambda_rest_bonds} {compound_thermodynamic_state_copy.lambda_alchemical_bonds_old} {compound_thermodynamic_state_copy.lambda_alchemical_bonds_new}")
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
        replica_id = np.where(self._replica_thermodynamic_states == 0)[0][0]
        print("ITERATION: ", self._iteration)
        print("REPLICA THERMOSTATES ", self._replica_thermodynamic_states, type(self._replica_thermodynamic_states))
        print("REPLICA ID ", replica_id, type(replica_id))
        self._reporter.write_sampler_states([self._sampler_states[replica_id]], self._iteration)
        
        self._reporter.write_replica_thermodynamic_states(self._replica_thermodynamic_states, self._iteration)
        self._reporter.write_mcmc_moves(self._mcmc_moves)  # MCMCMoves can store internal statistics.
        self._reporter.write_energies(self._energy_thermodynamic_states, self._neighborhoods, self._energy_unsampled_states,
                                      self._iteration)
        self._reporter.write_mixing_statistics(self._n_accepted_matrix, self._n_proposed_matrix, self._iteration)

# Set up sampler
_logger.setLevel(logging.DEBUG)
_logger.info("About to start repex")
_logger.info(f"move steps: {int((move_length*1000)/timestep)}")
_logger.info(f"timestep: {timestep} fs")
move = mcmc.LangevinSplittingDynamicsMove(timestep=timestep*unit.femtoseconds, n_steps=int((move_length*1000)/timestep), context_cache=context_cache)
simulation = ReplicaExchangeSampler2(mcmc_moves=move, number_of_iterations=length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{phase}_{name.lower()}_{length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()
