import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutorRBD
from simtk import unit

solvent_delivery = PointMutationExecutorRBD("0_rbd.pdb",
                        'R', # First chain is the barstar one
                        '439',
                        'LYS',
                        ligand_input="0_ace2.pdb",
                        ionic_strength=0.05*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True
                       )

outdir = "/data/chodera/zhangi/perses_benchmark/neq/14/test/"
apo_htf = solvent_delivery.get_apo_htf()
with open(os.path.join(outdir, "1_apo.pickle"), "wb") as f:
    pickle.dump(apo_htf, f)

complex_htf = solvent_delivery.get_complex_htf()
with open(os.path.join(outdir, "1_complex.pickle"), "wb") as f:
    pickle.dump(complex_htf, f)

apo_rhtf_0 = solvent_delivery.get_apo_rhtf_0()
with open(os.path.join(outdir, "1_apo_0.pickle"), "wb") as f:
    pickle.dump(apo_rhtf_0, f)

complex_rhtf_0 = solvent_delivery.get_complex_rhtf_0()
with open(os.path.join(outdir, "1_complex_0.pickle"), "wb") as f:
    pickle.dump(complex_rhtf_0, f)

apo_rhtf_1 = solvent_delivery.get_apo_rhtf_1()
with open(os.path.join(outdir, "1_apo_1.pickle"), "wb") as f:
    pickle.dump(apo_rhtf_1, f)

complex_rhtf_1 = solvent_delivery.get_complex_rhtf_1()
with open(os.path.join(outdir, "1_complex_1.pickle"), "wb") as f:
    pickle.dump(complex_rhtf_1, f)