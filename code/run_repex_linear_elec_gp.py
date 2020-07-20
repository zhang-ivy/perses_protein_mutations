import argparse
import pickle
import os
from perses.annihilation.lambda_protocol import LambdaProtocol
import simtk.unit as unit
from openmmtools.multistate import MultiStateReporter
from openmmtools import mcmc
import logging
import datetime

from perses.annihilation.lambda_protocol import RelativeAlchemicalState, LambdaProtocol

from openmmtools.multistate import sams, replicaexchange
from openmmtools import cache, utils
from perses.dispersed.utils import configure_platform
from openmmtools import cache
cache.global_context_cache.platform = configure_platform(utils.get_fastest_platform().getName())
from openmmtools.states import *
from perses.dispersed.utils import create_endstates
from openmmtools.integrators import LangevinIntegrator

import numpy as np
import copy
import logging

class HybridCompatibilityMixinWindows(object):
    """
    Mixin that allows the MultistateSampler to accommodate the situation where
    unsampled endpoints have a different number of degrees of freedom.
    """

    def __init__(self, *args, hybrid_factory=None, **kwargs):
        self._hybrid_factory = hybrid_factory
        super(HybridCompatibilityMixinWindows, self).__init__(*args, **kwargs)

    def setup(self, n_states, temperature, storage_file, minimisation_steps=100,
              n_replicas=None, lambda_schedule=None,
              lambda_protocol=LambdaProtocol(), endstates=True, unique_positions=None):


        from perses.dispersed import feptasks

        thermodynamic_state_list = []
        sampler_state_list = []

        context_cache = cache.ContextCache()

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
        
        for lambda_val, pos in zip(lambda_schedule, unique_positions):
            
            self._factory._new_positions = pos * unit.nanometer
            self._factory._hybrid_positions = self._factory._compute_hybrid_positions()
            positions = self._factory.hybrid_positions
            hybrid_system = self._factory.hybrid_system

            lambda_zero_alchemical_state = RelativeAlchemicalState.from_system(hybrid_system)

            thermostate = ThermodynamicState(hybrid_system, temperature=temperature)
            compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[lambda_zero_alchemical_state])
            sampler_state =  SamplerState(positions, box_vectors=hybrid_system.getDefaultPeriodicBoxVectors())
            
            compound_thermodynamic_state.set_alchemical_parameters(lambda_val,lambda_protocol)
            thermodynamic_state_list.append(compound_thermodynamic_state)

            # now generating a sampler_state for each thermodyanmic state, with relaxed positions
            integrator = LangevinIntegrator(temperature=temperature)
            context, context_integrator = context_cache.get_context(compound_thermodynamic_state, integrator)
            context.setPositions(positions)
            context.setPeriodicBoxVectors(*compound_thermodynamic_state.system.getDefaultPeriodicBoxVectors())
            feptasks.minimize(compound_thermodynamic_state, sampler_state)

            for i in range(5):
                try: 
                    context_integrator.step(100)
                    break
                except Exception as e:
                    print(e)
                    pass
            
            sampler_state_list.append(sampler_state)

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

class HybridRepexSamplerWindows(HybridCompatibilityMixinWindows, replicaexchange.ReplicaExchangeSampler):
    """
    ReplicaExchangeSampler that supports unsampled end states with a different number of positions
    """

    def __init__(self, *args, hybrid_factory=None, **kwargs):
        super(HybridRepexSamplerWindows, self).__init__(*args, hybrid_factory=hybrid_factory, **kwargs)
        self._factory = hybrid_factory


# Define lambda protocol
custom_functions = {'lambda_sterics_core':
                         lambda x: x,
                         'lambda_electrostatics_core':
                         lambda x: x,
                         'lambda_sterics_insert':
                         lambda x: 2.0 * x if x < 0.5 else 1.0,
                         'lambda_sterics_delete':
                         lambda x: 0.0 if x < 0.5 else 2.0 * (x - 0.5),
                         'lambda_electrostatics_insert':
                         lambda x: x,
                         'lambda_electrostatics_delete':
                         lambda x: x,
                         'lambda_bonds':
                         lambda x: x,
                         'lambda_angles':
                         lambda x: x,
                         'lambda_torsions':
                         lambda x: x
                         }


# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
args = parser.parse_args()

i = os.path.basename(os.path.dirname(args.dir))
apo_htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), "rb" ))

with open(os.path.join(args.dir, f"mmc2_barstar_T42_positions_complex.npy"), 'rb') as f:
    all_pos = np.load(f, allow_pickle=True)

# Build the hybrid repex samplers
suffix = 'run'; selection = 'not water'; checkpoint_interval = 10; n_states = 11; n_cycles = 5000
lambda_protocol = LambdaProtocol(functions=custom_functions)
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}.nc")
reporter = MultiStateReporter(reporter_file, analysis_particle_indices = apo_htf.hybrid_topology.select(selection), checkpoint_interval = checkpoint_interval)
hss = HybridRepexSamplerWindows(mcmc_moves=mcmc.LangevinSplittingDynamicsMove(timestep= 2.0 * unit.femtoseconds,
                                                                      collision_rate=5.0 / unit.picosecond,
                                                                      n_steps=500,
                                                                      reassign_velocities=False,
                                                                      n_restart_attempts=20,
                                                                      splitting="V R R R O R R R V",
                                                                      constraint_tolerance=1e-06),
                                                                      hybrid_factory=apo_htf, online_analysis_interval=10)
hss.setup(n_states=n_states, temperature=300*unit.kelvin, storage_file=reporter, lambda_protocol=lambda_protocol, endstates=False, unique_positions=all_pos[:11])
hss.extend(n_cycles)