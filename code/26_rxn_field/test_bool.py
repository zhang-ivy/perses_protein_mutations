import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from simtk import unit, openmm
import argparse
import os
import time
import mdtraj as md
from tqdm import tqdm
from openeye import oechem

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.INFO)

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('--is_rxn_field', dest='is_rxn_field', action='store_true', help='whether to use rxn field protocol')
parser.add_argument('--not_rxn_field', dest='is_rxn_field', action='store_false', help='whether to use rxn field protocol')
args = parser.parse_args()


# Define lambda functions
_logger.info(f"Defining lambda fn")

print(args.is_rxn_field)
x = 'lambda'
if not args.is_rxn_field:
    print("pme")

else:
    print("rxn field")
