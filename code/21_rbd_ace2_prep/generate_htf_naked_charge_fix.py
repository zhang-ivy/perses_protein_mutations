import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutorRBD
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping
from pkg_resources import resource_filename
from simtk import unit
import argparse

# Read args
parser = argparse.ArgumentParser(description='generate htf')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('residue', type=str, help='residue position')
parser.add_argument('mutant', type=str, help='three letter code for amino acid to mutate to')
parser.add_argument('--rbd_file', type=str, help='rbd file')
parser.add_argument('--ace2_file', type=str, help='ace2 file')
args = parser.parse_args()

rbd_file = resource_filename('perses', 'data/rbd-ace2/0_rbd.pdb') if not args.rbd_file else args.rbd_file
ace2_file = resource_filename('perses', 'data/rbd-ace2/0_ace2.pdb') if not args.ace2_file else args.ace2_file

import logging
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

from perses.utils.openeye import createOEMolFromSDF, extractPositionsFromOEMol, oechem
from perses.annihilation.relative import HybridTopologyFactory, RepartitionedHybridTopologyFactory
from perses.rjmc.topology_proposal import PointMutationEngine, PointMutationEngineRBD
from perses.rjmc.geometry import FFAllAngleGeometryEngine
from perses.utils.rbd import edit_pdb_for_tleap, edit_tleap_in_inputs, edit_tleap_in_ions, generate_tleap_system

import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as unit
import numpy as np
from openmoltools import forcefield_generators
import mdtraj as md
from openmmtools.constants import kB
from perses.tests.utils import validate_endstate_energies
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SystemGenerator
import os
from pkg_resources import resource_filename
import shutil
import tempfile

ENERGY_THRESHOLD = 1e-2
temperature = 298 * unit.kelvin
kT = kB * temperature
beta = 1.0/kT
ring_amino_acids = ['TYR', 'PHE', 'TRP', 'PRO', 'HIS']
CL_CHARGE = unit.Quantity(value=-1.0, unit=unit.elementary_charge)
CL_SIGMA = unit.Quantity(value=0.4477656957373345, unit=unit.nanometer)
CL_EPSILON = unit.Quantity(value=0.14891274399999999, unit=unit.kilojoule_per_mole)
NA_CHARGE = unit.Quantity(value=1.0, unit=unit.elementary_charge)
NA_SIGMA = unit.Quantity(value=0.2439280690268249, unit=unit.nanometer)
NA_EPSILON = unit.Quantity(value=0.3658460312, unit=unit.kilojoule_per_mole)
O_CHARGE = unit.Quantity(value=-0.834, unit=unit.elementary_charge)
H_CHARGE = unit.Quantity(value=0.417, unit=unit.elementary_charge)


class PointMutationExecutorRBD2(PointMutationExecutorRBD):
    def __init__(self,
                 protein_filename,
                 mutation_chain_id,
                 mutation_residue_id,
                 proposed_residue,
                 phase='complex',
                 clean=False,
                 conduct_endstate_validation=True,
                 ligand_input=None,
                 ligand_index=0,
                 water_model='tip3p',
                 ionic_strength=0.15 * unit.molar,
                 forcefield_files=['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml'],
                 barostat=openmm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 50),
                 forcefield_kwargs={'removeCMMotion': False, 'ewaldErrorTolerance': 0.00025, 'constraints' : app.HBonds, 'hydrogenMass' : 4 * unit.amus},
                 periodic_forcefield_kwargs={'nonbondedMethod': app.PME},
                 nonperiodic_forcefield_kwargs=None,
                 small_molecule_forcefields='gaff-2.11',
                 complex_box_dimensions=None,
                 apo_box_dimensions=None,
                 flatten_torsions=False,
                 flatten_exceptions=False,
                 vanilla=True,
                 repartitioned=True,
                 debug_dir=None,
                 **kwargs):
        """
        arguments
            protein_filename : str
                path to protein (to mutate); .pdb
            mutation_chain_id : str
                name of the chain to be mutated
            mutation_residue_id : str
                residue id to change
            proposed_residue : str
                three letter code of the residue to mutate to
            phase : str, default complex
                if phase == vacuum, then the complex will not be solvated with water; else, it will be solvated with tip3p
            clean : bool, default False
                whether to clean the PDB for tleap
            conduct_endstate_validation : bool, default True
                whether to conduct an endstate validation of the HybridTopologyFactory. If using the RepartitionedHybridTopologyFactory,
                endstate validation cannot and will not be conducted.
            ligand_file : str, default None
                path to ligand of interest (i.e. small molecule or protein); .sdf or .pdb
            ligand_index : int, default 0
                which ligand to use
            water_model : str, default 'tip3p'
                solvent model to use for solvation
            ionic_strength : float * unit.molar, default 0.15 * unit.molar
                the total concentration of ions (both positive and negative) to add using Modeller.
                This does not include ions that are added to neutralize the system.
                Note that only monovalent ions are currently supported.
            forcefield_files : list of str, default ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
                forcefield files for proteins and solvent
            barostat : openmm.MonteCarloBarostat, default openmm.MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin, 50)
                barostat to use
            forcefield_kwargs : dict, default {'removeCMMotion': False, 'ewaldErrorTolerance': 1e-4, 'constraints' : app.HBonds, 'hydrogenMass' : 4 * unit.amus}
                forcefield kwargs for system parametrization
            periodic_forcefield_kwargs : dict, default {'nonbondedMethod': app.PME}
                periodic forcefield kwargs for system parametrization
            nonperiodic_forcefield_kwargs : dict, default None
                non-periodic forcefield kwargs for system parametrization
            small_molecule_forcefields : str, default 'gaff-2.11'
                the forcefield string for small molecule parametrization
            complex_box_dimensions : Vec3, default None
                define box dimensions of complex phase;
                if None, padding is 1nm
            apo_box_dimensions :  Vec3, default None
                define box dimensions of apo phase phase;
                if None, padding is 1nm
            flatten_torsions : bool, default False
                in the htf, flatten torsions involving unique new atoms at lambda = 0 and unique old atoms are lambda = 1
            flatten_exceptions : bool, default False
                in the htf, flatten exceptions involving unique new atoms at lambda = 0 and unique old atoms at lambda = 1
            vanilla : bool, default True
                whether to generate a vanilla HybridTopologyFactory
            repartitioned : bool, default True
                whether to generate a RepartitionedHybridTopologyFactory
            debug_dir : str, default None
                if specified, debug output files will be saved here
        TODO : allow argument for spectator ligands besides the 'ligand_file'
        """
        # Make debug directory
        is_temp = False
        if debug_dir:
            if not os.path.exists(debug_dir):
                os.system(f"mkdir {debug_dir}")
        else:
            debug_dir = tempfile.mkdtemp()
            is_temp = True        

        ## Generate the old topology, positions, and system
        # Prep PDBs for tleap
        _logger.info("Editing PDBs for tleap")
        protein_name = os.path.basename(protein_filename)
        ligand_name = os.path.basename(ligand_input)
        protein_tleap = os.path.join(debug_dir, f"{protein_name[:-4]}_tleap.pdb")
        ligand_tleap = os.path.join(debug_dir, f"{ligand_name[:-4]}_tleap.pdb")
        if clean:
            edit_pdb_for_tleap(protein_filename, protein_tleap)
            edit_pdb_for_tleap(ligand_input, ligand_tleap)
        else:
            os.system(f"cp {protein_filename} {protein_tleap}")
            os.system(f"cp {ligand_input} {ligand_tleap}")
        
        # Edit tleap files
        _logger.info("Editing tleap.in input files")
        apo_tleap_prefix = os.path.join(debug_dir, "1_rbd_tleap")
        complex_tleap_prefix = os.path.join(debug_dir, "1_rbd_ace2_tleap")
        apo_template = resource_filename('perses', 'data/rbd-ace2/1_rbd_template_tleap.in')
        complex_template = resource_filename('perses', 'data/rbd-ace2/1_rbd_ace2_template_tleap.in')
        edit_tleap_in_inputs(apo_template, apo_tleap_prefix, debug_dir)
        edit_tleap_in_inputs(complex_template, complex_tleap_prefix, debug_dir)
 
        _logger.info("Editing tleap.in number of ions")
        edit_tleap_in_ions(apo_tleap_prefix)
        edit_tleap_in_ions(complex_tleap_prefix)
        
        # Generating old systems
        _logger.info("Generating solvated old systems")
        apo_topology, apo_positions, apo_system = generate_tleap_system(apo_tleap_prefix)
        complex_topology, complex_positions, complex_system = generate_tleap_system(complex_tleap_prefix)
        
        # Correct the topologies
        _logger.info("Correcting tleap topologies")
        apo_topology_corrected = self._correct_topology(apo_topology)
        complex_topology_corrected = self._correct_topology(complex_topology, is_apo=False)
        
        # Format inputs for pipeline
        inputs = [[apo_topology_corrected, apo_positions, apo_system, apo_tleap_prefix, False], [complex_topology_corrected, complex_positions, complex_system, complex_tleap_prefix, True]]
    
        # Make system generator -- note this is only for system_generator.forcefield call in PointMutationEngine init
        molecules = []
        self.system_generator = SystemGenerator(forcefields=forcefield_files,
                                                barostat=barostat,
                                                forcefield_kwargs=forcefield_kwargs,
                                                periodic_forcefield_kwargs=periodic_forcefield_kwargs,
                                                nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs,
                                                small_molecule_forcefield=small_molecule_forcefields,
                                                molecules=molecules,
                                                cache=None)
        
        # Run pipeline...
        htfs = []
        for (top, pos, sys, tleap_prefix, is_complex) in inputs:
            name = 'rbd_ace2' if is_complex else 'rbd'
            _logger.info(f"Generating topology proposal for {name}")
            point_mutation_engine = PointMutationEngineRBD(wildtype_topology=top,
                                                         system_generator=self.system_generator,
                                                         chain_id=mutation_chain_id, # Denote the chain id allowed to mutate (it's always a string variable)
                                                         max_point_mutants=1,
                                                         residues_allowed_to_mutate=[mutation_residue_id], # The residue ids allowed to mutate
                                                         allowed_mutations=[(mutation_residue_id, proposed_residue)], # The residue ids allowed to mutate with the three-letter code allowed to change
                                                         aggregate=True) # Always allow aggregation

            topology_proposal, new_positions = point_mutation_engine.propose(sys, top, pos, tleap_prefix, is_complex, debug_dir)

            # Fix naked charges in old and new systems
            for system in [topology_proposal.old_system, topology_proposal.new_system]:
                force_dict = {i.__class__.__name__: i for i in system.getForces()}
                if 'NonbondedForce' in [i for i in force_dict.keys()]:
                    nb_force = force_dict['NonbondedForce']
                    for i in range(nb_force.getNumParticles()):
                        charge, sigma, epsilon = nb_force.getParticleParameters(i)
                        if sigma == 0*unit.nanometer:
                            sigma = 0.6*unit.nanometer
                            nb_force.setParticleParameters(i, charge, sigma, epsilon)
                        if epsilon == 0*unit.kilojoule_per_mole:
                            epsilon = 0.0001*unit.kilojoule_per_mole
                            nb_force.setParticleParameters(i, charge, sigma, epsilon)
            
            # Check for charge change...
            charge_diff = point_mutation_engine._get_charge_difference(current_resname = topology_proposal._old_topology.residue_topology.name,
                                                                       new_resname = topology_proposal._new_topology.residue_topology.name)
            _logger.info(f"charge diff: {charge_diff}")
            if charge_diff != 0:
                new_water_indices_to_ionize = point_mutation_engine.get_water_indices(charge_diff, new_positions, topology_proposal._new_topology, radius=0.8)
                _logger.info(f"new water indices to ionize {new_water_indices_to_ionize}")
                PointMutationExecutorRBD._modify_new_system(new_water_indices_to_ionize, topology_proposal._new_system, charge_diff)
                PointMutationExecutorRBD._modify_atom_classes(new_water_indices_to_ionize, topology_proposal)

            factories = []
            if vanilla:
                repartitioned_endstate = None
                self.generate_htf(HybridTopologyFactory, topology_proposal, pos, new_positions, flatten_exceptions, flatten_torsions, repartitioned_endstate, is_complex)
            if repartitioned:
                for repartitioned_endstate in [0, 1]:
                    self.generate_htf(RepartitionedHybridTopologyFactory, topology_proposal, pos, new_positions, flatten_exceptions, flatten_torsions, repartitioned_endstate, is_complex)


solvent_delivery = PointMutationExecutorRBD2(rbd_file,
                        'R',
                        args.residue,
                        args.mutant,
                        ligand_input=ace2_file,
                        ionic_strength=0.05*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True, 
                        debug_dir=os.path.join(args.dir, "debug/")
                       )

outdir = args.dir
i = os.path.basename(os.path.dirname(args.dir))

apo_htf = solvent_delivery.get_apo_htf()
with open(os.path.join(outdir, f"{i}_apo.pickle"), "wb") as f:
    pickle.dump(apo_htf, f)

complex_htf = solvent_delivery.get_complex_htf()
with open(os.path.join(outdir, f"{i}_complex.pickle"), "wb") as f:
    pickle.dump(complex_htf, f)

apo_rhtf_0 = solvent_delivery.get_apo_rhtf_0()
with open(os.path.join(outdir, f"{i}_apo_0.pickle"), "wb") as f:
    pickle.dump(apo_rhtf_0, f)

complex_rhtf_0 = solvent_delivery.get_complex_rhtf_0()
with open(os.path.join(outdir, f"{i}_complex_0.pickle"), "wb") as f:
    pickle.dump(complex_rhtf_0, f)

apo_rhtf_1 = solvent_delivery.get_apo_rhtf_1()
with open(os.path.join(outdir, f"{i}_apo_1.pickle"), "wb") as f:
    pickle.dump(apo_rhtf_1, f)

complex_rhtf_1 = solvent_delivery.get_complex_rhtf_1()
with open(os.path.join(outdir, f"{i}_complex_1.pickle"), "wb") as f:
    pickle.dump(complex_rhtf_1, f)

# Render atom map
atom_map_filename = f'{outdir}/atom_map.png'
render_protein_residue_atom_mapping(apo_htf._topology_proposal, atom_map_filename)