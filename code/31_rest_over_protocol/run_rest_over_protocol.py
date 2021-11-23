import pickle
import numpy as np
import argparse
import os
import logging

from simtk import openmm
from simtk.openmm import unit

import openmmtools
from openmmtools.states import SamplerState, ThermodynamicState, CompoundThermodynamicState
from openmmtools.multistate import MultiStateReporter
from openmmtools import mcmc
from openmmtools import cache
from openmmtools.constants import kB

from perses.annihilation.lambda_protocol import RESTState, RESTCapableRelativeAlchemicalState
from perses.dispersed import feptasks

parser = argparse.ArgumentParser(description='run rest over protocol repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='phase')
parser.add_argument('n_states', type=int, help='num of states')
parser.add_argument('n_cycles', type=int, help='num of iterations to run')
parser.add_argument('T_max', type=int, help='T_max for rest')
args = parser.parse_args()

# Load pickled htf
i = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), "rb" ))
hybrid_system = htf.hybrid_system
hybrid_positions = htf.hybrid_positions
box_vectors = hybrid_system.getDefaultPeriodicBoxVectors()

# Subclass the HybridCompatilityMixin and HybridRepexSampler classes

#from perses.annihilation.lambda_protocol import RelativeAlchemicalState, LambdaProtocol
from openmmtools.multistate import sams, replicaexchange
from openmmtools import cache, utils
from perses.dispersed.utils import configure_platform
cache.global_context_cache.platform = configure_platform(utils.get_fastest_platform().getName())
context_cache = cache.ContextCache(capacity=None, time_to_live=None)
from openmmtools.states import CompoundThermodynamicState, SamplerState, ThermodynamicState
from perses.dispersed.utils import create_endstates

import numpy as np
import copy

class HybridCompatibilityMixin(object):
    """
    Mixin that allows the MultistateSampler to accommodate the situation where
    unsampled endpoints have a different number of degrees of freedom.
    """

    def __init__(self, *args, hybrid_factory=None, **kwargs):
        self._hybrid_factory = hybrid_factory
        super(HybridCompatibilityMixin, self).__init__(*args, **kwargs)

    def setup(self, n_states, temperature, T_max, storage_file, minimisation_steps=100,
              n_replicas=None, lambda_schedule=None, endstates=True):


        from perses.dispersed import feptasks

        hybrid_system = self._factory.hybrid_system
        positions = self._factory.hybrid_positions

        # lambda_zero_alchemical_state = RelativeAlchemicalState.from_system(hybrid_system)
        # thermostate = ThermodynamicState(hybrid_system, temperature=temperature)
        # compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[lambda_zero_alchemical_state])

        # Create thermodynamic state
        lambda_zero_alchemical_state = RESTCapableRelativeAlchemicalState.from_system(hybrid_system)
        thermostate = ThermodynamicState(hybrid_system, temperature=temperature)
        compound_thermodynamic_state = CompoundThermodynamicState(thermostate,
                                                                  composable_states=[lambda_zero_alchemical_state])

        # Set alchemical parameters
        beta_0 = 1 / (kB * temperature)
        beta_m = 1 / (kB * T_max)
        global_lambda = 0
        compound_thermodynamic_state.set_alchemical_parameters(global_lambda, beta_0, beta_m)

        thermodynamic_state_list = []
        sampler_state_list = []

#         context_cache = cache.ContextCache()

        if n_replicas is None:
            _logger.info(f'n_replicas not defined, setting to match n_states, {n_states}')
            n_replicas = n_states
        elif n_replicas > n_states:
            _logger.warning(f'More sampler states: {n_replicas} requested greater than number of states: {n_states}. Setting n_replicas to n_states: {n_states}')
            n_replicas = n_states

        # TODO this feels like it should be somewhere else... just not sure where. Maybe into lambda_protocol
        if lambda_schedule is None:
            lambda_schedule = np.linspace(0.,1.,n_states)
        else:
            assert (len(lambda_schedule) == n_states) , 'length of lambda_schedule must match the number of states, n_states'
            assert (lambda_schedule[0] == 0.), 'lambda_schedule must start at 0.'
            assert (lambda_schedule[-1] == 1.), 'lambda_schedule must end at 1.'
            difference = np.diff(lambda_schedule)
            assert ( all(i >= 0. for i in difference ) ), 'lambda_schedule must be monotonicly increasing'

        #starting with the initial positions generated py geometry.py
        sampler_state =  SamplerState(positions, box_vectors=hybrid_system.getDefaultPeriodicBoxVectors())
        for lambda_val in lambda_schedule:
            compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
            compound_thermodynamic_state_copy.set_alchemical_parameters(lambda_val, beta_0, beta_m)
            # compound_thermodynamic_state_copy.set_alchemical_parameters(lambda_val,lambda_protocol)
            thermodynamic_state_list.append(compound_thermodynamic_state_copy)

             # now generating a sampler_state for each thermodyanmic state, with relaxed positions
            context, context_integrator = context_cache.get_context(compound_thermodynamic_state_copy)
            feptasks.minimize(compound_thermodynamic_state_copy,sampler_state, max_iterations=0)
            sampler_state_list.append(copy.deepcopy(sampler_state))

        reporter = storage_file

        # making sure number of sampler states equals n_replicas
        if len(sampler_state_list) != n_replicas:
            # picking roughly evenly spaced sampler states
            # if n_replicas == 1, then it will pick the first in the list
            idx = np.round(np.linspace(0, len(sampler_state_list) - 1, n_replicas)).astype(int)
            sampler_state_list = [state for i,state in enumerate(sampler_state_list) if i in idx]

        assert len(sampler_state_list) == n_replicas

        if endstates:
            # generating unsampled endstates
            _logger.info('Generating unsampled endstates.')
            unsampled_dispersion_endstates = create_endstates(copy.deepcopy(thermodynamic_state_list[0]), copy.deepcopy(thermodynamic_state_list[-1]))
            self.create(thermodynamic_states=thermodynamic_state_list, sampler_states=sampler_state_list,
                    storage=reporter, unsampled_thermodynamic_states=unsampled_dispersion_endstates)
        else:
            self.create(thermodynamic_states=thermodynamic_state_list, sampler_states=sampler_state_list,
                        storage=reporter)


        def print_global_parameters(state):
            print("REST")
            print(state.lambda_rest_bonds)
            print(state.lambda_rest_angles)
            print(state.lambda_rest_torsions)
            print(state.lambda_rest_electrostatics)
            print(state.lambda_rest_electrostatics_exceptions)
            print(state.lambda_rest_sterics)
            print(state.lambda_rest_sterics_exceptions)
            print("ALCHEMICAL")
            print(state.lambda_alchemical_bonds_old)
            print(state.lambda_alchemical_bonds_new)
            print(state.lambda_alchemical_angles_old)
            print(state.lambda_alchemical_angles_new)
            print(state.lambda_alchemical_torsions_old)
            print(state.lambda_alchemical_torsions_new)
            print(state.lambda_alchemical_electrostatics_old)
            print(state.lambda_alchemical_electrostatics_new)
            print(state.lambda_alchemical_electrostatics_exceptions_old)
            print(state.lambda_alchemical_electrostatics_exceptions_new)
            print(state.lambda_alchemical_sterics_old)
            print(state.lambda_alchemical_sterics_new)
            print(state.lambda_alchemical_sterics_exceptions_old)
            print(state.lambda_alchemical_sterics_exceptions_new)
            print(state.lambda_alchemical_electrostatics_reciprocal)

        for i, state in enumerate(thermodynamic_state_list):
            print(i)
            print_global_parameters(state)

class HybridRepexSampler(HybridCompatibilityMixin, replicaexchange.ReplicaExchangeSampler):
    """
    ReplicaExchangeSampler that supports unsampled end states with a different number of positions
    """

    def __init__(self, *args, hybrid_factory=None, **kwargs):
        super(HybridRepexSampler, self).__init__(*args, hybrid_factory=hybrid_factory, **kwargs)
        self._factory = hybrid_factory

# Instantiate sampler 
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
reporter_file = os.path.join(os.path.join(args.dir, f"{i}_{args.phase}.nc"))
reporter = MultiStateReporter(reporter_file, checkpoint_interval=10)
hss = HybridRepexSampler(mcmc_moves=mcmc.LangevinSplittingDynamicsMove(timestep= 4.0 * unit.femtoseconds,
                                                                      collision_rate=1.0 / unit.picosecond,
                                                                      n_steps=250,
                                                                      reassign_velocities=False,
                                                                      n_restart_attempts=20,
                                                                      splitting="V R R R O R R R V",
                                                                      constraint_tolerance=1e-06, 
                                                                      context_cache=context_cache),
                                                                      replica_mixing_scheme='swap-all',
                                                                      hybrid_factory=htf, 
                                                                      online_analysis_interval=10)
hss.setup(n_states=args.n_states, temperature=300*unit.kelvin, T_max=args.T_max * unit.kelvin, storage_file=reporter, endstates=False)

# Run simulation
hss.extend(args.n_cycles)
