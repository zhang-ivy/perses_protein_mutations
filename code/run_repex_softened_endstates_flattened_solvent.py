from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState
from openmmtools import states
from simtk import openmm, unit
from openmmtools.mcmc import LangevinSplittingDynamicsMove, GHMCMove
from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import numpy as np
import pickle
import argparse
import os

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA, corresponding to the endstate')
parser.add_argument('endstate', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('direction', type=str, help="forward or reverse")
args = parser.parse_args()

if args.direction == "forward":
	# Generate htf for capped ALA->THR in vacuum
	atp, sys_gen = generate_atp(phase="solvent")

	# At alanine endstate
	htf = generate_dipeptide_top_pos_sys(atp.topology, 
	                                         new_res='THR', 
	                                         system=atp.system, 
	                                         positions=atp.positions,
	                                         system_generator=sys_gen, 
	                                         conduct_htf_prop=True,
	                                         repartitioned=True,
	                                         endstate=args.endstate,
	                                         flatten_torsions=True,
	                                         flatten_exceptions=True,
	                                         validate_endstate_energy=False)

elif args.direction == 'reverse':
	# Generate htf for capped THR->ALA in vacuum
	pdb = app.PDBFile("../input/thr_vacuum.pdb")

	forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
	barostat = openmm.MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin, 50)
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
	modeller = app.Modeller(atp.topology, atp.positions)
	modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=9*unit.angstroms, ionicStrength=0.15*unit.molar)
	solvated_topology = modeller.getTopology()
	solvated_positions = modeller.getPositions()

	# Canonicalize the solvated positions: turn tuples into np.array
	positions = unit.quantity.Quantity(value=np.array([list(atom_pos) for atom_pos in solvated_positions.value_in_unit_system(unit.md_unit_system)]), unit=unit.nanometers)
	system = system_generator.create_system(solvated_topology)

	htf = generate_dipeptide_top_pos_sys(solvated_topology, 
	                                         new_res='ALA', 
	                                         system=system, 
	                                         positions=positions,
	                                         system_generator=system_generator, 
	                                         conduct_htf_prop=True, 
	                                         repartitioned=True,
	                                         endstate=args.endstate,
	                                         flatten_torsions=True,
	                                         flatten_exceptions=True,
	                                         validate_endstate_energy=False)

# Render atom map
atom_map_filename = f'{args.dir}/atom_map.png'
render_protein_residue_atom_mapping(htf._topology_proposal, atom_map_filename)

# Save htf as pickle
i = os.path.basename(os.path.dirname(args.dir))
pickle.dump(htf, open(os.path.join(args.dir, f"{i}_{args.phase}_{args.name}.pickle"), "wb" ))

# Alchemify the hybrid system
atoms_to_alchemify = list(htf._atom_classes['unique_new_atoms']) + list(htf._atom_classes['unique_old_atoms'])

alch_factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
alchemical_region = AlchemicalRegion(alchemical_atoms=list(atoms_to_alchemify), alchemical_torsions=True)
alchemical_system = alch_factory.create_alchemical_system(htf.hybrid_system, alchemical_region)

# Initialize compound thermodynamic states at different temperatures and alchemical states.
protocol = {'temperature': [300]*unit.kelvin*11,
            'lambda_electrostatics': [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0],
            'lambda_sterics': [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0],
           'lambda_torsions': [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]}
alchemical_state = AlchemicalState.from_system(alchemical_system)
compound_states = states.create_thermodynamic_state_protocol(alchemical_system, 
                                                             protocol=protocol, 
                                                             composable_states=[alchemical_state])

# Set up the sampler
n_steps = 500 # 1 ps
n_iterations = args.length*1000 

# # Propagate the replicas with Langevin dynamics.
# langevin_move = LangevinSplittingDynamicsMove(timestep=2.0*unit.femtosecond, n_steps=n_steps)
# simulation = ReplicaExchangeSampler(mcmc_moves=langevin_move, number_of_iterations=n_iterations)

# Propagate the replicas with GHMC move.
ghmc_move = GHMCMove(timestep=2.0*unit.femtosecond, n_steps=n_steps)
simulation = ReplicaExchangeSampler(mcmc_moves=ghmc_move, number_of_iterations=n_iterations)

# Run simulation
i = os.path.basename(os.path.dirname(args.dir))
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=states.SamplerState(htf.hybrid_positions),
                  storage=reporter)
simulation.run()
