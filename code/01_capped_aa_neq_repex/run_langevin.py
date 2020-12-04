import logging
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
import pickle
import numpy as np
import time
from openmmtools.integrators import LangevinIntegrator
from simtk import unit
from simtk import openmm

with open("/data/chodera/zhangi/perses_benchmark/repex/0/solvent/ALA_CYS_solvent.pickle", 'rb') as f:
    htf = pickle.load(f)

# Set up Langevin integrator
nsteps_eq = 62500 # 0.25 ns
nsteps_neq= 20000 # 80 ps
neq_splitting='V R H O R V'
timestep=4.0 * unit.femtosecond
integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, timestep)

# Set up context
platform_name = 'CUDA'
system = htf.hybrid_system
positions = htf.hybrid_positions
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')

context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run 500 steps of MD
md_step_times = []
nsteps = 500
for step in range(nsteps):
    initial_time = time.time()
    integrator.step(step)
    elapsed_time = (time.time() - initial_time) * unit.seconds
    md_step_times.append(elapsed_time)
    _logger.info(f'step: {step} took: {elapsed_time / unit.seconds} seconds')

md_step_times = [time / unit.seconds for time in md_step_times]

from matplotlib import pyplot as plt
fig = plt.figure(dpi=150)
fig.patch.set_facecolor('white')
plt.plot(md_step_times)
fig.savefig("langevin_step_times.png")

pickle.dump(md_step_times, open("langevin_step_times.pickle", "wb" ))

