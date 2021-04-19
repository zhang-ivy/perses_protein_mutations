import math
from simtk import unit
from openmmtools import testsystems, states, mcmc, multistate
from openmmtools.multistate import ReplicaExchangeSampler
testsystem = testsystems.AlanineDipeptideExplicit()
import os
import tempfile
import logging

n_replicas = 11  # Number of temperature replicas.
T_min = 300.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0) for i in range(n_replicas)]
thermodynamic_states = [states.ThermodynamicState(system=testsystem.system, temperature=T) for T in temperatures]

move = mcmc.LangevinSplittingDynamicsMove(timestep=4.0*unit.femtoseconds, n_steps=250, collision_rate=5.0 / unit.picosecond, reassign_velocities=False, n_restart_attempts=20, constraint_tolerance=1e-06)
simulation = ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=5000, replica_mixing_scheme='swap-neighbors', online_analysis_interval=10)

storage_path = '/data/chodera/zhangi/perses_benchmark/repex/31/0/3/alanine_dipeptide.nc'
reporter = multistate.MultiStateReporter(storage_path, checkpoint_interval=10)
simulation.create(thermodynamic_states=thermodynamic_states, sampler_states=states.SamplerState(testsystem.positions, box_vectors=testsystem.system.getDefaultPeriodicBoxVectors()), storage=reporter)

_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

simulation.run()

