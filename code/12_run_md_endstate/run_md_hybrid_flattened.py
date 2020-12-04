import logging
import pickle
import numpy as np
from openmmtools.integrators import LangevinIntegrator
from openmmtools.states import ThermodynamicState, CompoundThermodynamicState
from perses.annihilation.lambda_protocol import RelativeAlchemicalState, LambdaProtocol
from simtk import openmm, unit
import argparse
import os
import time
from simtk.openmm.app import PDBFile
import mdtraj as md
from tqdm import tqdm
from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('endstate', type=int, help="0 or 1")
parser.add_argument('direction', type=str, help="forward or reverse")
args = parser.parse_args()

# Define simulation parameters
temperature = 300 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
nsteps = 10000000 # 20 ns
timestep = 2.0 * unit.femtosecond
platform_name = 'CUDA'

if args.direction == "forward":
    # Generate htf for capped ALA->THR in vacuum
    atp, sys_gen = generate_atp()

    # At alanine endstate
    htf = generate_dipeptide_top_pos_sys(atp.topology, 
                                             new_res='THR', 
                                             system=atp.system, 
                                             positions=atp.positions,
                                             system_generator=sys_gen, 
                                             conduct_htf_prop=True,
                                             flatten_torsions=True,
                                             flatten_exceptions=True,
                                             validate_endstate_energy=False)

elif args.direction == 'reverse':
    # Generate htf for capped THR->ALA in vacuum
    pdb = app.PDBFile("../input/thr_vacuum.pdb")

    forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    barostat = None
    system_generator = SystemGenerator(forcefields=forcefield_files,
                                   barostat=barostat,
                                   forcefield_kwargs={'removeCMMotion': False,
                                                        'ewaldErrorTolerance': 1e-4,
                                                        'constraints' : app.HBonds,
                                                        'hydrogenMass' : 4 * unit.amus},
                                    periodic_forcefield_kwargs=None,
                                    small_molecule_forcefield='gaff-2.11',
                                    nonperiodic_forcefield_kwargs = {'nonbondedMethod': app.NoCutoff}, 
                                       molecules=None, 
                                       cache=None)
    system = system_generator.create_system(pdb.topology) 
    positions = unit.quantity.Quantity(value = np.array([list(atom_pos) for atom_pos in pdb.positions.value_in_unit_system(unit.md_unit_system)]), unit = unit.nanometers)

    htf = generate_dipeptide_top_pos_sys(pdb.topology, 
                                             new_res='ALA', 
                                             system=system, 
                                             positions=positions,
                                             system_generator=system_generator, 
                                             conduct_htf_prop=True,
                                             flatten_torsions=True,
                                             flatten_exceptions=True,
                                             validate_endstate_energy=False)

# Render atom map
atom_map_filename = f'{args.dir}/atom_map.png'
render_protein_residue_atom_mapping(htf._topology_proposal, atom_map_filename)

i = os.path.basename(os.path.dirname(args.dir))
pickle.dump(htf, open(os.path.join(args.dir, f"{i}_{args.phase}_.pickle"), "wb" ))

system = htf.hybrid_system
positions = htf.hybrid_positions

# Create CompoundThermodynamicState at the appropriate endstate
lambda_alchemical_state = RelativeAlchemicalState.from_system(system)
lambda_protocol = LambdaProtocol(functions = 'default')
lambda_alchemical_state.set_alchemical_parameters(args.endstate, lambda_protocol)
thermodynamic_state = CompoundThermodynamicState(ThermodynamicState(system, temperature=temperature), composable_states=[lambda_alchemical_state])

# Set up integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = thermodynamic_state.create_context(integrator, platform=platform)
context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)

# Minimize
openmm.LocalEnergyMinimizer.minimize(context)

# Run equilibration
final_pos = np.empty(shape=(40001, htf.hybrid_topology.n_atoms, 3))
pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
i = 0
final_pos[i] = pos * unit.nanometers
for step in tqdm(range(nsteps)):
    initial_time = time.time()
    integrator.step(1)
    if step % 2500 == 0:
        pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
        final_pos[i] = pos *unit.nanometers
        i += 1
        elapsed_time = (time.time() - initial_time) * unit.seconds
        _logger.info(f'Step: {step} took {elapsed_time} seconds')

# Save traj
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}_equil_hybrid_{args.endstate}_20ns.npy"), 'wb') as f:
    np.save(f, final_pos)
