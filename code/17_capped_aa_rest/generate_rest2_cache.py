import pickle
import os
from perses.annihilation.rest import RESTTopologyFactory
from perses.annihilation.lambda_protocol import RESTStateV2
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

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('T_max', type=int, help='in kelvin')
args = parser.parse_args()

# Load rhtf
i = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}_{args.state}.pickle"), "rb" ))

# Build REST factory
if args.phase == 'solvent':
	# factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(6, 20)) + [0, 1, 20, 22])
	factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(6, 20)) + [1549, 1550, 1551, 1552]) # THR + off ALA atoms
elif args.phase == 'vacuum':
	factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(6, 20)) + [26, 27, 28, 29]) # THR + off ALA atoms


# Get REST system
REST_system = factory.REST_system

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = args.T_max * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

# Create reference thermodynamic state
lambda_zero_alchemical_state = RESTStateV2.from_system(REST_system)
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
    feptasks.minimize(compound_thermodynamic_state_copy, sampler_state)
    sampler_state_list.append(copy.deepcopy(sampler_state))

# Set up sampler
move = mcmc.GHMCMove(timestep=4.0*unit.femtoseconds, n_steps=250)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()