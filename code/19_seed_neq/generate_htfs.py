import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutor
import argparse
from simtk import openmm, unit

# Set up logger
import logging
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run t-repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('resid', type=str, help='resid of residue to be mutated')
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('sim_number', type=int, help='index of job array, starts at 1')
parser.add_argument('input_file', type=str, default="../../input/mmc2_barstar.pdb", help='name of input file ')
parser.add_argument('ligand_file', type=str, default="../../input/mmc2_barnase.pdb", help='name of ligand file')
args = parser.parse_args()

if args.sim_number == 1:
    phase = 'apo'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 2:
    phase = 'apo'
    name = args.new_aa_name
    state = 1
elif args.sim_number == 3:
    phase = 'complex'
    name = args.old_aa_name
    state = 0
elif args.sim_number == 4:
    phase = 'complex'
    name = args.new_aa_name
    state = 1

# Load rhtf or generate it if it doesn't exist
i = os.path.basename(os.path.dirname(args.dir))
path = os.path.join(args.dir, f"{i}_{phase}_{state}.pickle")
if not os.path.exists(path):
    solvent_delivery = PointMutationExecutor(args.input_file,
                        '1', # First chain is the barstar one
                        args.resid,
                        args.new_aa_name,
                        ligand_file=args.ligand_file,
                        ionic_strength=0.05*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True,
                        conduct_endstate_validation=False
                       )
    pickle.dump(solvent_delivery.get_apo_htf(), open(os.path.join(args.dir, f"{i}_apo_{state}.pickle"), "wb" ))
    pickle.dump(solvent_delivery.get_complex_htf(), open(os.path.join(args.dir, f"{i}_complex_{state}.pickle"), "wb" ))