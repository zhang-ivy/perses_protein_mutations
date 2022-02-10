import logging
import pickle
import numpy as np
from openmmtools.integrators import PeriodicNonequilibriumIntegrator
from openmmtools.constants import kB
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
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='solvent or vacuum')
parser.add_argument('sim_number', type=int, help='number in job name - 1')
parser.add_argument('old_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('new_aa_name', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('length', type=float, help='neq switching time in ns')
parser.add_argument('endstate', type=int, help='starting endstate (0 or 1)')
parser.add_argument('T_max', type=int, help='t_max for rest, in Kelvin')
parser.add_argument('chain_A', type=int, help='(first) chain index for which to to add a virtual bond for complex phase')
parser.add_argument('chain_B', type=int, help='(second) chain index for which to to add a virtual bond for complex phase')
parser.add_argument('protocol', type=str, help='protocol to use: default or bridge')
parser.add_argument('--cache', type=int, default=1, help='length of rest cache in ns')
parser.add_argument('--restrained', action='store_true', help='whether to use restrained cache')
args = parser.parse_args()

# Define simulation parameters
nsteps_eq = 10
nsteps_neq = int(args.length*250000) # 1 ns
neq_splitting='V R H O R V'
timestep = 4.0 * unit.femtosecond
platform_name = 'CUDA'
cache_length = args.cache if args.cache else 1
temperature = 300.0 * unit.kelvin
aa_name = args.old_aa_name if args.endstate == 0 else args.new_aa_name

# Define lambda functions
x = '(1 - lambda)'
beta0 = (1 / (kB * temperature)).value_in_unit_system(unit.md_unit_system)
beta = (1 / (kB * args.T_max * unit.kelvin)).value_in_unit_system(unit.md_unit_system)

if args.protocol == 'bridge':
    ALCHEMICAL_FUNCTIONS = {
                          'lambda_rest_bonds': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_angles': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_torsions': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_electrostatics': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_electrostatics_exceptions': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_sterics': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_sterics_exceptions': f"select(step({x} - 0.25), select(step({x} - 0.75), 4 * (1 - sqrt({beta} / {beta0})) * ({x} - 1) + 1, sqrt({beta}/{beta0})), -4 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
#'lambda_rest_bonds': "1",
#'lambda_rest_angles': "1",
#'lambda_rest_torsions':"1",
#'lambda_rest_electrostatics': "1",
#'lambda_rest_electrostatics_exceptions': "1",
#'lambda_rest_sterics':"1",
#'lambda_rest_sterics_exceptions': "1", 

                         'lambda_alchemical_bonds_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_bonds_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_angles_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_angles_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_torsions_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_torsions_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_electrostatics_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_electrostatics_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_electrostatics_exceptions_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_electrostatics_exceptions_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_electrostatics_reciprocal': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_sterics_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_sterics_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         'lambda_alchemical_sterics_exceptions_old': f'select(step({x} - 0.25), select(step({x} - 0.75), 0.0, 1.5 - 2 * {x}), 1.0)',
                         'lambda_alchemical_sterics_exceptions_new': f'select(step({x} - 0.25), select(step({x} - 0.75), 1.0, 2 * {x} - 0.5), 0.0)',
                         }

elif args.protocol == 'default':
    ALCHEMICAL_FUNCTIONS = {
                          'lambda_rest_bonds': f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_angles': f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_torsions':f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_electrostatics': f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_electrostatics_exceptions': f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_sterics':f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
                         'lambda_rest_sterics_exceptions': f"select(step({x} - 0.5), 2 * (1 - sqrt({beta} / {beta0})) * {x} - 1 + 2 * sqrt({beta} / {beta0}), -2 * (1 - sqrt({beta} / {beta0})) * {x} + 1)",
#'lambda_rest_bonds': "1",
#'lambda_rest_angles': "1",
#'lambda_rest_torsions':"1",
#'lambda_rest_electrostatics': "1",
#'lambda_rest_electrostatics_exceptions': "1",
#'lambda_rest_sterics':"1",
#'lambda_rest_sterics_exceptions': "1", 

                         'lambda_alchemical_bonds_old': f'1 - {x}',
                         'lambda_alchemical_bonds_new': x,
                         'lambda_alchemical_angles_old': f'1 - {x}',
                         'lambda_alchemical_angles_new': x,
                         'lambda_alchemical_torsions_old': f'1 - {x}',
                         'lambda_alchemical_torsions_new': x,
                         'lambda_alchemical_electrostatics_old': f'1 - {x}',
                         'lambda_alchemical_electrostatics_new': x,
                         'lambda_alchemical_electrostatics_exceptions_old': f'1 - {x}',
                         'lambda_alchemical_electrostatics_exceptions_new': x,
                         'lambda_alchemical_electrostatics_reciprocal': x,
                         'lambda_alchemical_sterics_old': f'1 - {x}',
                         'lambda_alchemical_sterics_new': x,
                         'lambda_alchemical_sterics_exceptions_old': f'1 - {x}',
                         'lambda_alchemical_sterics_exceptions_new': x
                         }
else:
    raise Exception("protocol must be 'default' or 'bridge'")


# Read in vanilla htf
i = os.path.basename(os.path.dirname(args.dir))
with open(os.path.join(args.dir, f"{i}_{args.phase}.pickle"), 'rb') as f:
    htf = pickle.load(f)
system = htf.hybrid_system

# Make sure LRC is set correctly
system.getForce(5).setUseLongRangeCorrection(False)
system.getForce(8).setUseDispersionCorrection(True)
_logger.info(f"CustomNonbondedForce_sterics use LRC? {system.getForce(5).getUseLongRangeCorrection()}")
_logger.info(f"NonbondedForce_sterics use LRC? {system.getForce(8).getUseDispersionCorrection()}")

# Add virtual bond for complex phase
if args.phase == 'complex':
    chains = list(htf.hybrid_topology.chains)
    atom_A = list(chains[args.chain_A].atoms)[0]
    atom_B = list(chains[args.chain_B].atoms)[0]
    force = openmm.CustomBondForce('0')
    force.addBond(atom_A.index, atom_B.index, [])
    system.addForce(force)
    _logger.info(f"Added virtual bond between {atom_A} and {atom_B}")

# Read in cache
append = '_restrained' if args.restrained else ''
_logger.info(f"Using restrained cache? {args.restrained}")
with open(os.path.join(args.dir, f"{i}_{args.phase}_{aa_name}_{cache_length}ns_snapshots{append}.npy"), 'rb') as f:
    cache_pos = np.load(f)
positions = cache_pos[args.sim_number]

# Read in cache box vectors
with open(os.path.join(args.dir, f"{i}_{args.phase}_{aa_name}_{cache_length}ns_box_vectors{append}.npy"), 'rb') as f:
    cache_box_vectors = np.load(f)
box_vectors = cache_box_vectors[args.sim_number][0]

# Set up integrator
integrator = PeriodicNonequilibriumIntegrator(ALCHEMICAL_FUNCTIONS, nsteps_eq, nsteps_neq, neq_splitting, timestep=timestep, temperature=temperature)

# Set up context
platform = openmm.Platform.getPlatformByName(platform_name)
if platform_name in ['CUDA', 'OpenCL']:
    platform.setPropertyDefaultValue('Precision', 'mixed')
if platform_name in ['CUDA']:
    platform.setPropertyDefaultValue('DeterministicForces', 'true')
context = openmm.Context(system, integrator, platform)
context.setPeriodicBoxVectors(*box_vectors)
context.setPositions(positions)
context.setVelocitiesToTemperature(temperature)
integrator.reset()

#for k, v in context.getParameters().items():
#    if 'old' in k:
#        context.setParameter(k, 0.0)
#    elif 'new' in k or 'reciprocal' in k:
#        context.setParameter(k, 1.0)

_logger.info("before eq")
for k, v in context.getParameters().items():
    _logger.info(f"{k} {v}")

# Run eq
integrator.step(nsteps_eq)
_logger.info("finished eq")

for k, v in context.getParameters().items():
    _logger.info(f"{k} {v}")

# Run neq
works_master = list()
neq_old, neq_new = list(), list()
neq_old_waters, neq_new_waters = list(), list()
works = [integrator.get_protocol_work(dimensionless=True)]
for step in range(int(nsteps_neq / 2500)):
    integrator.step(2500)
    _logger.info(f"Neq: {step*2500} completed")
    
#    for k, v in context.getParameters().items():
#        _logger.info(f"{k} {v}")
#    for i in range(integrator.getNumGlobalVariables()):
#        name = integrator.getGlobalVariableName(i)
#        _logger.info(f"{name}: {integrator.getGlobalVariable(i)}")
#    _logger.info("\n")

    works.append(integrator.get_protocol_work(dimensionless=True))
    
    pos = context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)
    old_pos = np.asarray(htf.old_positions(pos))
    old_traj = md.Trajectory(old_pos, md.Topology.from_openmm(htf._topology_proposal.old_topology))
    old_pos_solute = old_traj.atom_slice(old_traj.top.select("not water")).xyz[0]
    
    new_pos = np.asarray(htf.new_positions(pos))
    new_traj = md.Trajectory(new_pos, md.Topology.from_openmm(htf._topology_proposal.new_topology))
    new_pos_solute = new_traj.atom_slice(new_traj.top.select("not water")).xyz[0]
    
    neq_old.append(old_pos_solute)
    neq_new.append(new_pos_solute)
    if step == int(nsteps_neq / 2500) - 1:
        neq_old_waters.append(old_pos)
        neq_new_waters.append(new_pos)
works_master.append(works)

for k, v in context.getParameters().items():
    _logger.info(f"{k} {v}")

# Save works
i = os.path.basename(os.path.dirname(args.dir))
direction = 'forward' if args.endstate == 0 else 'reverse'
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_{direction}.npy"), 'wb') as f:
    np.save(f, works_master)

with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_{direction}_neq_old.npy"), 'wb') as f:
    np.save(f, neq_old)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_{direction}_neq_new.npy"), 'wb') as f:
    np.save(f, neq_new)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_{direction}_neq_old_waters.npy"), 'wb') as f:
    np.save(f, neq_old_waters)
with open(os.path.join(args.dir, f"{i}_{args.phase}_{args.sim_number}_{direction}_neq_new_waters.npy"), 'wb') as f:
    np.save(f, neq_new_waters)
