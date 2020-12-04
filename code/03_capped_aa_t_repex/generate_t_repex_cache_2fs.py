import argparse
import math
from simtk import unit
from openmmtools import testsystems, states, mcmc, multistate
import os
import tempfile
import pickle
from perses.annihilation.lambda_protocol import RelativeAlchemicalState
import mdtraj as md
import numpy as np
from simtk.unit.quantity import Quantity
from perses.tests.test_topology_proposal import generate_atp, generate_dipeptide_top_pos_sys
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping


# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('state', type=int, help='aka lambda, e.g. 0 or 1')
parser.add_argument('length', type=int, help='in ns')
parser.add_argument('direction', type=str, help="forward or reverse")
args = parser.parse_args()

i = os.path.basename(os.path.dirname(args.dir))

if args.direction == "forward":
	# Generate htf for capped ALA->THR in vacuum
	atp, sys_gen = generate_atp()

	# At alanine endstate
	htf = generate_dipeptide_top_pos_sys(atp.topology, 
	                                         new_res = 'THR', 
	                                         system = atp.system, 
	                                         positions = atp.positions,
	                                         system_generator = sys_gen, 
	                                         conduct_htf_prop=True,
	                                         flatten_torsions=True, 
	                                         validate_endstate_energy=False)
	pickle.dump(htf, open(os.path.join(args.dir, f"{i}_vacuum.pickle"), "wb" ))

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
	                                         new_res = 'ALA', 
	                                         system = system, 
	                                         positions = positions,
	                                         system_generator = system_generator, 
	                                         conduct_htf_prop=True, 
	                                         flatten_torsions=True, 
	                                         validate_endstate_energy=False)
	pickle.dump(htf, open(os.path.join(args.dir, f"{i}_vacuum.pickle"), "wb" ))

# Render atom map
atom_map_filename = f'{args.dir}/atom_map.png'
render_protein_residue_atom_mapping(htf._topology_proposal, atom_map_filename)

# Create states for each replica
n_replicas = 12  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 600.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

alchemical_state = RelativeAlchemicalState.from_system(htf.hybrid_system)
alchemical_state.set_alchemical_parameters(args.state)

thermodynamic_states = [states.ThermodynamicState(system=htf.hybrid_system, temperature=T)
                        for T in temperatures]
compound_states = [states.CompoundThermodynamicState(thermodynamic_state=thermodynamic_state, composable_states=[alchemical_state])
                        for thermodynamic_state in thermodynamic_states]

# Set up sampler
move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=500)
simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=args.length*1000)

# Run t-repex
reporter_file = os.path.join(args.dir, f"{i}_{args.phase}_{args.name.lower()}_{args.length}ns.nc")
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=1)
simulation.create(thermodynamic_states=compound_states,
                  sampler_states=states.SamplerState(htf.hybrid_positions),
                  storage=reporter)
simulation.run()