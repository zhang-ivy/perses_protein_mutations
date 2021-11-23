import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutor
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping
from simtk import unit
import mdtraj as md
import argparse

# Read args
parser = argparse.ArgumentParser(description='generate htf')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('residue', type=str, help='residue position')
parser.add_argument('mutant', type=str, help='three letter code for amino acid to mutate to')
parser.add_argument('--rbd_file', type=str, help='rbd file')
parser.add_argument('--ace2_file', type=str, help='ace2 file')
args = parser.parse_args()

rbd_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/rbd_protonated.pdb"  if not args.rbd_file else args.rbd_file
ace2_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/ace2_protonated.pdb" if not args.ace2_file else args.ace2_file


import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as unit
import numpy as np

# Set up logger
import logging
_logger = logging.getLogger("setup")
_logger.setLevel(logging.INFO)

class PointMutationExecutor2(PointMutationExecutor):
    def _solvate(self,
               topology,
               positions,
               water_model,
               phase,
               ionic_strength,
               box_dimensions=None):
        """
        Generate a solvated topology, positions, and system for a given input topology and positions.
        For generating the system, the forcefield files provided in the constructor will be used.
        Parameters
        ----------
        topology : app.Topology
            Topology of the system to solvate
        positions : [n, 3] ndarray of Quantity nm
            the positions of the unsolvated system
        forcefield : SystemGenerator.forcefield
            forcefield file of solvent to add
        water_model : str
            solvent model to use for solvation
        phase : str
            if phase == vacuum, then the complex will not be solvated with water; else, it will be solvated with tip3p
        ionic_strength : float * unit.molar
            the total concentration of ions (both positive and negative) to add using Modeller.
            This does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        Returns
        -------
        solvated_topology : app.Topology
            Topology of the system with added waters
        solvated_positions : [n + 3(n_waters), 3] ndarray of Quantity nm
            Solvated positions
        solvated_system : openmm.System
            The parameterized system, containing a barostat if one was specified.
        """
        modeller = app.Modeller(topology, positions)
    
        geompadding = 0.9 * unit.nanometers
        maxSize = max(max((pos[i] for pos in positions))-min((pos[i] for pos in positions)) for i in range(3))
        vectors = openmm.Vec3(1,0,0), openmm.Vec3(1/3,2*np.sqrt(2)/3,0), openmm.Vec3(-1/3,np.sqrt(2)/3,np.sqrt(6)/3)
        boxVectors = [(maxSize+geompadding)*v for v in vectors]
        
        # Now we have to add missing atoms
        if phase != 'vacuum':
            _logger.info(f"solvating at {ionic_strength} using {water_model}")
            modeller.addSolvent(self.system_generator.forcefield, model=water_model, boxVectors=boxVectors, ionicStrength=ionic_strength)
        else:
            pass

        solvated_topology = modeller.getTopology()
        solvated_positions = modeller.getPositions()

        # Canonicalize the solvated positions: turn tuples into np.array
        solvated_positions = unit.quantity.Quantity(value=np.array([list(atom_pos) for atom_pos in solvated_positions.value_in_unit_system(unit.md_unit_system)]), unit=unit.nanometers)
        solvated_system = self.system_generator.create_system(solvated_topology)

        return solvated_topology, solvated_positions, solvated_system


solvent_delivery = PointMutationExecutor2(rbd_file,
                        '1',
                        args.residue,
                        args.mutant,
                        ligand_input=ace2_file,
                        forcefield_files=['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', '/home/zhangi/choderalab/og_openmmforcefields/openmmforcefields/amber/ffxml/GLYCAM_06j-1.xml'],
                        ionic_strength=0.15*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True,
                        generate_unmodified_hybrid_topology_factory=True,
                        generate_rest_capable_hybrid_topology_factory=True,
                        conduct_endstate_validation=False
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

# Save pdbs
traj = md.Trajectory(complex_htf.old_positions(complex_htf.hybrid_positions), md.Topology.from_openmm(complex_htf._topology_proposal.old_topology))
traj.save(os.path.join(outdir, "complex_old.pdb"))

traj = md.Trajectory(complex_htf.new_positions(complex_htf.hybrid_positions), md.Topology.from_openmm(complex_htf._topology_proposal.new_topology))
traj.save(os.path.join(outdir, "complex_new.pdb"))

traj = md.Trajectory(apo_htf.old_positions(apo_htf.hybrid_positions), md.Topology.from_openmm(apo_htf._topology_proposal.old_topology))
traj.save(os.path.join(outdir, "apo_old.pdb"))

traj = md.Trajectory(apo_htf.new_positions(apo_htf.hybrid_positions), md.Topology.from_openmm(apo_htf._topology_proposal.new_topology))
traj.save(os.path.join(outdir, "apo_new.pdb"))
