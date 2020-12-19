import pickle
import os
from perses.annihilation.rest import RESTTopologyFactory
from perses.annihilation.lambda_protocol import RESTState
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
import numpy as np
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import pickle

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('T_max', type=int, help='in kelvin')
args = parser.parse_args()

# Make test system
if args.name == 'ALA':
    # from openmmtools.testsystems import AlanineDipeptideExplicit
    # ala = AlanineDipeptideExplicit()
    # sys = ala.system
    # sys.removeForce(4)
    # positions = ala.positions
    pdb = app.PDBFile("../../input/ala_vacuum.pdb")
elif args.name == 'THR':
    pdb = app.PDBFile("../../input/thr_vacuum.pdb")

forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
barostat = openmm.MonteCarloBarostat(1.0 * unit.atmosphere, 298 * unit.kelvin, 50)
system_generator = SystemGenerator(forcefields=forcefield_files,
                               barostat=barostat,
                               forcefield_kwargs={'removeCMMotion': False,
                                                    'ewaldErrorTolerance': 1e-4,
                                                    'constraints' : app.HBonds,
                                                    'hydrogenMass' : 4 * unit.amus},
                                periodic_forcefield_kwargs={'nonbondedMethod': app.PME},
                                small_molecule_forcefield='gaff-2.11',
                                nonperiodic_forcefield_kwargs=None, 
                                   molecules=None, 
                                   cache=None)
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=9*unit.angstroms, ionicStrength=0.15*unit.molar)
solvated_topology = modeller.getTopology()
solvated_positions = modeller.getPositions()

# Canonicalize the solvated positions: turn tuples into np.array
positions = unit.quantity.Quantity(value=np.array([list(atom_pos) for atom_pos in solvated_positions.value_in_unit_system(unit.md_unit_system)]), unit=unit.nanometers)
sys = system_generator.create_system(solvated_topology)

# Save vanilla topology and system
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.name.lower()}_vanilla_system.pickle"), "wb") as f:
    pickle.dump(sys, f)
with open(os.path.join(args.dir, f"{i}_{args.name.lower()}_vanilla_topology.pickle"), "wb") as f:
    pickle.dump(solvated_topology, f)

# Build REST factory
factory = RESTTopologyFactory(sys, solute_region=list(range(22)))

# Get REST system
REST_system = factory.REST_system

# Create states for each replica
temperatures = [300, 420, 600]

# Create reference thermodynamic state
lambda_zero_alchemical_state = RESTState.from_system(REST_system)
thermostate = ThermodynamicState(REST_system, temperature=T_min)
compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[lambda_zero_alchemical_state])

context_cache = cache.ContextCache()

# Create thermodynamics states
sampler_state =  SamplerState(positions, box_vectors=sys.getDefaultPeriodicBoxVectors())
beta_0 = 1/(kB*T_min)
thermodynamic_state_list = []
sampler_state_list = []
for temperature in temperatures:
    beta_m = 1/(kB*temperature)
    compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
    compound_thermodynamic_state_copy.set_alchemical_parameters(beta_0, beta_m)
    thermodynamic_state_list.append(compound_thermodynamic_state_copy)

    # now generating a sampler_state for each thermodyanmic state, with relaxed positions
    context, context_integrator = context_cache.get_context(compound_thermodynamic_state_copy)
    feptasks.minimize(compound_thermodynamic_state_copy, sampler_state)
    sampler_state_list.append(copy.deepcopy(sampler_state))

# Set up sampler
move = mcmc.GHMCMove(timestep=4.0*unit.femtoseconds, n_steps=250)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000,  replica_mixing_scheme='swap-neighbors')

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()