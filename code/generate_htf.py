import argparse
import itertools
import pickle
from perses.app.relative_point_mutation_setup import PointMutationExecutor
import simtk.openmm.app as app
import os

# Read WT and proposed residue abbreviations
parser = argparse.ArgumentParser(description='generate htfs')
parser.add_argument('output_dir', type=str, help='path to output htf')
args = parser.parse_args()


amino_acids = ['ALA', 'CYS', 'SER', 'THR']
# Solvent
for wt, proposed in itertools.permutations(amino_acids,r=2):
     # Create hybrid topology factory
     apm_delivery = PointMutationExecutor(f"../input/{wt.lower()}_vacuum.pdb", 
                              '1', 
                              '2', 
                              proposed,
                             )
     apo_htf = apm_delivery.get_apo_htf()
     pickle.dump(apo_htf, open(os.path.join(args.output_dir, f"solvent/{wt}_{proposed}_solvent.pickle"), "wb" ))

# Vacuum
for wt, proposed in itertools.permutations(amino_acids,r=2):
     # Create hybrid topology factory
     apm_delivery = PointMutationExecutor(f"../input/{wt.lower()}_vacuum.pdb", 
                              '1', 
                              '2', 
                              proposed,
                              phase='vacuum',
	                          barostat=None,
	                          periodic_forcefield_kwargs=None, 
	                          nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
                             )
     apo_htf = apm_delivery.get_apo_htf()
     pickle.dump(apo_htf, open(os.path.join(args.output_dir, f"vacuum/{wt}_{proposed}_vacuum.pickle"), "wb" ))
