import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutor
import argparse
from simtk import openmm, unit
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping

# Set up logger
import logging
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('resid', type=str, help='resid of residue to be mutated')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('--input_file', type=str, default="/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/1brs_barnase_renumbered.pdb", help='name of input file ')
parser.add_argument('--ligand_file', type=str, default="/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/1brs_barstar_renumbered.pdb", help='name of ligand file')
args = parser.parse_args()

_logger.info(f"input file: {args.input_file}")
_logger.info(f"ligand file: {args.ligand_file}")

# Load rhtf or generate it if it doesn't exist
solvent_delivery = PointMutationExecutor(args.input_file,
                        '1', # First chain is the barstar one
                        args.resid,
                        args.new_aa_name,
                        ligand_input=args.ligand_file,
                        ionic_strength=0.05*unit.molar,
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
