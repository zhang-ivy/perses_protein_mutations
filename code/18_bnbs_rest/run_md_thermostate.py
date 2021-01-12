import pickle
import os
from perses.annihilation.rest import RESTTopologyFactory
from perses.annihilation.lambda_protocol import RESTStateV2
from openmmtools.states import SamplerState, ThermodynamicState, CompoundThermodynamicState
from simtk import openmm, unit
import math
from openmmtools.constants import kB
from openmmtools import mcmc, multistate
import copy
from perses.dispersed import feptasks
import mdtraj as md
import numpy as np
from openmmtools.integrators import LangevinIntegrator
from tqdm import tqdm

# Load in r-htf at endstate = 0
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/12/0"
htf = pickle.load(open(os.path.join(out_dir, f"0_apo_0.pickle"), "rb" ))

# Build REST factory
query_indices = list(range(669, 683)) + list(range(14877, 14881))
traj = md.Trajectory(np.array(htf.hybrid_positions), htf.hybrid_topology)
rest_atoms = list(md.compute_neighbors(traj, 0.3, query_indices)[0])
# water_atoms = list(md.compute_neighbors(traj, 0.8, query_indices, haystack_indices=list(range(1441, htf.hybrid_topology.n_atoms)))[0])
factory = RESTTopologyFactory(htf.hybrid_system, solute_region=list(set(rest_atoms)))

# Get REST system
REST_system = factory.REST_system

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 1200.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

# Create reference thermodynamic state
lambda_zero_alchemical_state = RESTStateV2.from_system(REST_system)
thermostate = ThermodynamicState(REST_system, temperature=T_min)
compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[lambda_zero_alchemical_state])

# Create thermodynamics states
sampler_state =  SamplerState(htf.hybrid_positions, box_vectors=htf.hybrid_system.getDefaultPeriodicBoxVectors())
beta_0 = 1/(kB*T_min)
thermodynamic_state_list = []
sampler_state_list = []
for temperature in temperatures:
    beta_m = 1/(kB*temperature)
    compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
    compound_thermodynamic_state_copy.set_alchemical_parameters(beta_0, beta_m)
#     print("solute", compound_thermodynamic_state_copy.solute_scale)
#     print("inter", compound_thermodynamic_state_copy.inter_scale)
    thermodynamic_state_list.append(compound_thermodynamic_state_copy)
    
    # now generating a sampler_state for each thermodynamic state, with relaxed positions
    # context, context_integrator = context_cache.get_context(compound_thermodynamic_state_copy)
    feptasks.minimize(compound_thermodynamic_state_copy, sampler_state)
    sampler_state_list.append(copy.deepcopy(sampler_state))

# Define simulation parameters
temperature = 298 * unit.kelvin
nsteps = 25000 # 100 ps
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
collision_rate = 1.0 / unit.picoseconds

# Set up integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')

context = thermodynamic_state_list[11].create_context(integrator, platform=platform)    

context.setPeriodicBoxVectors(*htf.hybrid_system.getDefaultPeriodicBoxVectors())
context.setPositions(htf.hybrid_positions)
context.setVelocitiesToTemperature(temperature)
print(compound_thermodynamic_state.reduced_potential(context))

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run equilibration
pos_old = []
for step in tqdm(range(nsteps)):
    integrator.step(1)
    if step % 25 == 0:
	    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
	    old_pos = np.asarray(htf.old_positions(pos))
	    pos_old.append(old_pos)

traj_old = md.Trajectory(pos_old, md.Topology.from_openmm(htf._topology_proposal.old_topology))
traj_old.save(os.path.join(out_dir, f"state_11_md_old.dcd")) 

traj_old[0].save(os.path.join(out_dir, f"state_11_md_old.pdb"))
