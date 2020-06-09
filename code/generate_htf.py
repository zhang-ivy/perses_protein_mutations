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
for i, (wt, proposed) in enumerate(itertools.permutations(amino_acids,r=2)):
    
    print(f"solvent: {wt}, {proposed}")
    if not os.path.exists(os.path.join(args.output_dir, f"{i}/")):
        os.makedirs(os.path.join(args.output_dir, f"{i}/"))

    # Create hybrid topology factory
    apm_delivery = PointMutationExecutor(f"../input/{wt.lower()}_vacuum.pdb", 
                            '1', 
                            '2', 
                            proposed,
                           )
    apo_htf = apm_delivery.get_apo_htf()
    pickle.dump(apo_htf, open(os.path.join(args.output_dir, f"{i}/{i}_solvent.pickle"), "wb" ))


# Vacuum
for i, (wt, proposed) in enumerate(itertools.permutations(amino_acids,r=2)):
    print(f"vacuum: {wt}, {proposed}")
    if not os.path.exists(os.path.join(args.output_dir, f"{i}/")):
        os.makedirs(os.path.join(args.output_dir, f"{i}/"))

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
    pickle.dump(apo_htf, open(os.path.join(args.output_dir, f"{i}/{i}_vacuum.pickle"), "wb" ))


