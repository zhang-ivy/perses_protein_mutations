import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutorRBD
from simtk import unit
import argparse

# Read args
parser = argparse.ArgumentParser(description='generate htf')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('residue', type=str, help='apo or complex')
parser.add_argument('mutant', type=str, help='number in job name - 1')
args = parser.parse_args()

solvent_delivery = PointMutationExecutorRBD("0_rbd.pdb",
                        'R',
                        args.residue,
                        args.mutant,
                        ligand_input="0_ace2.pdb",
                        ionic_strength=0.05*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True
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