from simtk import openmm, unit
from openmmtools.integrators import PeriodicNonequilibriumIntegrator, LangevinIntegrator, NonequilibriumLangevinIntegrator, AlchemicalNonequilibriumLangevinIntegrator
from perses.dispersed.utils import configure_platform
from openmmtools import utils, cache, states
import os
import pickle

# Load htf
dir_num = 5
with open(f"5_solvent.pickle", "rb") as f:
    apo_htf = pickle.load(f)

system = apo_htf.hybrid_system
positions = apo_htf.hybrid_positions
box_vectors = apo_htf.hybrid_system.getDefaultPeriodicBoxVectors()

# Define simulation params
nsteps_eq = 25000 # 100 ps
nsteps_neq = 25000 # 100 ps
neq_splitting='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
temperature = 300.0 * unit.kelvin

# Define lambda protocol
x = 'lambda'
ALCHEMICAL_FUNCTIONS = {'lambda_0_bonds_old': f'1 - {x}',
                         'lambda_0_bonds_new': x,
                         'lambda_0_angles_old': f'1 - {x}',
                         'lambda_0_angles_new': x,
                         'lambda_0_torsions_old': f'1 - {x}',
                         'lambda_0_torsions_new': x,
                         'lambda_0_electrostatics_old': f'1 - {x}',
                         'lambda_0_electrostatics_new': x,
                         'lambda_0_electrostatics_exceptions_old': f'1 - {x}', 
                         'lambda_0_electrostatics_exceptions_new': x,
                         'lambda_0_sterics_old': f'1 - {x}',
                         'lambda_0_sterics_new': x,
                         'lambda_0_sterics_exceptions_old': f'1 - {x}',
                         'lambda_0_sterics_exceptions_new': x}

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep, temperature=temperature)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*box_vectors)
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run eq forward (0 -> 1)
for step in range(int(nsteps_eq/10)):
    integrator.step(10)
    print(step)


# Run neq forward (0 -> 1)
for fwd_step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    print(f"Forward neq: {fwd_step*2500} completed")

# Run eq reverse (1 -> 0)
for step in range(nsteps_eq):
    integrator.step(1)
    print(step)

