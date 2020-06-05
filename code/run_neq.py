import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from simtk import unit
from simtk import openmm

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Define lambda functions
x = 'lambda'
DEFAULT_ALCHEMICAL_FUNCTIONS = {
                             'lambda_sterics_core': x,
                             'lambda_electrostatics_core': x,
                             'lambda_sterics_insert': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_sterics_delete': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_insert': f"select(step({x} - 0.5), 2.0 * ({x} - 0.5), 0.0)",
                             'lambda_electrostatics_delete': f"select(step({x} - 0.5), 1.0, 2.0 * {x})",
                             'lambda_bonds': x,
                             'lambda_angles': x,
                             'lambda_torsions': x}

# Define simulation parameters
nsteps_eq = 62500 # 0.25 ns
nsteps_neq= 20000 # 80 ps
neq_splitting='V R H O R V'
timestep=4.0 * unit.femtosecond
temperature = 300*unit.kelvin
collision_rate = 90/unit.picosecond
platform_name = 'CUDA'
system = htf.hybrid_system
positions = htf.hybrid_positions

# Read in htf
with open("/data/chodera/zhangi/perses_benchmark/repex/0/solvent/ALA_CYS_solvent.pickle", 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system
positions = htf.hybrid_positions

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(DEFAULT_ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPositions(positions)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run neq
ncycles = 100
forward_works, reverse_works = list(), list()
for _ in range(ncycles):
    # Equilibrium (lambda = 0)
    integrator.step(nsteps_eq)
    
    # Forward (0 -> 1)
    forward_works.append(integrator.get_protocol_work(dimensionless=True))
    for fwd_step in range(nsteps_neq):
        integrator.step(fwd_step)
        forward_works.append(integrator.get_protocol_work(dimensionless=True))
        
    # Equilibrium (lambda = 1)
    integrator.step(nsteps_eq)
    
    # Reverse work (1 -> 0)
    reverse_works.append(integrator.get_protocol_work(dimensionless=True))
    for rev_step in range(nsteps_neq):
        integrator.step(rev_step)
        reverse_works.append(integrator.get_protocol_work(dimensionless=True))
        
# Save works
np.save("ALA_CYS_solvent_forward.npy", forward_works)
np.save("ALA_CYS_solvent_reverse.npy", reverse_works)