import argparse
import itertools
import pickle
from perses.app.relative_point_mutation_setup import PointMutationExecutor
import simtk.openmm.app as app
import os

amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                        'SER', 'THR', 'TRP', 'TYR', 'VAL']
outdir = "/data/chodera/zhangi/perses_benchmark/test_maps/"

# Vacuum
for wt in amino_acids[1:]:
    if wt != 'ALA':
    for proposed in amino_acids:
        if wt != proposed:
            print(f"vacuum: {wt}, {proposed}")
            # Create hybrid topology factory
            apo_delivery = PointMutationExecutor(os.path.join(outdir, f"{wt.lower()}_vacuum.pdb"), 
                                      '1', 
                                      '2', 
                                      proposed,
                                      phase='vacuum',
                                  barostat=None,
                                  periodic_forcefield_kwargs=None, 
                                  nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
                                     )
            pickle.dump(apo_delivery.get_apo_htf(), open(os.path.join(outdir, f"{wt.lower()}_{proposed.lower()}_vacuum.pickle"), "wb" ))


