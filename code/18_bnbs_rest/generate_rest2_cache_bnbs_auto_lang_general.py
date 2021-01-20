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
parser.add_argument('input_file', type=str, default="../../input/mmc2_barstar.pdb", help='name of input file ')
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
radius = 0.2

# Load rhtf or generate it if it doesn't exist
i = os.path.basename(os.path.dirname(args.dir))
path = os.path.join(args.dir, f"{i}_{phase}_{state}.pickle")
if not os.path.exists(path):
    solvent_delivery = PointMutationExecutor(args.input_file,
                        '1', # First chain is the barstar one
                        args.resid,
                        args.new_aa_name,
                        ligand_file="../../input/mmc2_barnase.pdb",
                        ionic_strength=0.05*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True,
                        repartitioned_endstate=state
                       )
    pickle.dump(solvent_delivery.get_apo_htf(), open(os.path.join(args.dir, f"{i}_apo_{state}.pickle"), "wb" ))
    pickle.dump(solvent_delivery.get_complex_htf(), open(os.path.join(args.dir, f"{i}_complex_{state}.pickle"), "wb" ))
htf = pickle.load(open(path, "rb" ))

# Build REST factory
_logger.info("Generating REST factory")
_logger.info(f"radius:{radius} nm")
query_indices = [atom.index for atom in list(htf.hybrid_topology.residues)[int(args.resid)].atoms]
_logger.info(f"query indices {query_indices}")
traj = md.Trajectory(np.array(htf.hybrid_positions), htf.hybrid_topology)
solute_atoms = list(traj.topology.select("is_protein"))
rest_atoms = list(md.compute_neighbors(traj, radius, query_indices, haystack_indices=solute_atoms)[0])
factory = RESTTopologyFactory(htf.hybrid_system, solute_region=rest_atoms)

_logger.info("Generating REST states")
# Get REST system
REST_system = factory.REST_system

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 1200 * unit.kelvin  # Maximum temperature.
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
_logger.setLevel(logging.DEBUG)
_logger.info("About to start repex")
print(f"move steps: {int((move_length*1000)/timestep)}")
print(f"timestep: {timestep} fs")
move = mcmc.LangevinSplittingDynamicsMove(timestep=timestep*unit.femtoseconds, n_steps=int((move_length*1000)/timestep))
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{phase}_{name.lower()}_{length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()