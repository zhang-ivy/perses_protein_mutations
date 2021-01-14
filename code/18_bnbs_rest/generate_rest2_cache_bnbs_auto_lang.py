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
import mdtraj as md
import numpy as np

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('T_max', type=int, help='in kelvin')
parser.add_argument('direction', type=str, help='forward or backward', default='forward')
parser.add_argument('move_length', type=float, help='length (in ps) of LangevinSplittingDynamicsMove')
parser.add_argument('timestep', type=float, help='2 or 4 fs')
parser.add_argument('radius', type=float, help='e.g. cutoff in nm for REST2 region')
parser.add_argument('states', type=int, help='for repex')
args = parser.parse_args()

# Load rhtf
i = os.path.basename(os.path.dirname(args.dir))
htf = pickle.load(open(os.path.join(args.dir, f"{i}_{args.phase}_{args.state}.pickle"), "rb" ))

# Build REST factory
# if args.phase == 'apo' and args.direction == 'backward':
#     # traj = md.Trajectory(np.array(htf.hybrid_positions), htf.hybrid_topology)
#     # rest_atoms = md.compute_neighbors(traj, 0.5, list(range(669, 683)))[0]
#     factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(669, 683)) + list(range(14877, 14881)))
# elif args.phase == 'complex' and args.direction == 'backward':
#     factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(669, 683)) + list(range(29474, 29478)))
# elif args.phase == 'apo' and args.direction == 'forward':
#     factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(669, 679)) + list(range(14873, 14881)))
# elif args.phase == 'complex' and args.direction == 'forward':
#     factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(range(669, 679)) + list(range(29470, 29478)))

if args.phase == 'apo' and args.direction == 'backward':
    query_indices = list(range(669, 683)) + list(range(14877, 14881))
elif args.phase == 'complex' and args.direction == 'backward':
    query_indices = list(range(669, 683)) + list(range(29474, 29478))
elif args.phase == 'apo' and args.direction == 'forward':
    query_indices = list(range(669, 679)) + list(range(14873, 14881))
elif args.phase == 'complex' and args.direction == 'forward':
    query_indices = list(range(669, 679)) + list(range(29470, 29478))

print(f"radius:{args.radius}")
traj = md.Trajectory(np.array(htf.hybrid_positions), htf.hybrid_topology)
solute_atoms = list(traj.topology.select("is_protein"))
rest_atoms = list(md.compute_neighbors(traj, args.radius, query_indices, haystack_indices=solute_atoms)[0])
# water_atoms = list(md.compute_neighbors(traj, 0.8, query_indices, haystack_indices=list(range(1441, htf.hybrid_topology.n_atoms)))[0])
factory = RESTTopologyFactory(htf.hybrid_system, solute_region=rest_atoms)

# Get REST system
REST_system = factory.REST_system

# Create states for each replica
n_replicas = args.states  # Number of temperature replicas.
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
print(f"move steps: {int((args.move_length*1000)/args.timestep)}")
print(f"timestep: {args.timestep}")
move = mcmc.LangevinSplittingDynamicsMove(timestep=args.timestep*unit.femtoseconds, n_steps=int((args.move_length*1000)/args.timestep))
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()