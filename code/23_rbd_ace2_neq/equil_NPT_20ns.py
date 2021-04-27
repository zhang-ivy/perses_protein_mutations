from __future__ import division, print_function

import sys
from sys import stdout

import mdtraj as md
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from mdtraj.reporters import DCDReporter
from simtk.openmm import XmlSerializer
from simtk.openmm import MonteCarloBarostat
from simtk.openmm.app import CheckpointReporter, PDBFile
from simtk.openmm.app.amberinpcrdfile import AmberInpcrdFile
from simtk.openmm.app.amberprmtopfile import AmberPrmtopFile
from simtk.openmm.app.pdbreporter import PDBReporter
from simtk.openmm.app.statedatareporter import StateDataReporter
import openmmtools
import time
from simtk.openmm import CustomBondForce
import argparse
import os

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('prefix', type=str, help='prefix to tleap')
args = parser.parse_args()

# Set parameters
print("Reading parameters...")
pressure = 1.0 * unit.atmospheres
temperature = 310 * unit.kelvin
nonbonded_method = app.PME
constraints = app.HBonds
remove_cm_motion = False
hydrogen_mass = 4.0 * unit.amu  # Using HMR
collision_rate = 1.0 / unit.picoseconds
timestep = 0.004 * unit.picoseconds  # We can use a 4fs timestep with HMR

# Set steps and frequencies
nsteps = 5000000  # 20 ns
report_freq = 100
chk_freq = 500
traj_freq = 4000  # 1250 frames

# Set input file names
amber_prmtop_file = f"{args.prefix}_tleap.prmtop"
amber_inpcrd_file = f"{args.prefix}_tleap.inpcrd"

# Set file names
output_prefix = f'/data/chodera/zhangi/perses_benchmark/neq/15/75/equil/{os.path.basename(args.prefix)}/'
integrator_xml_filename = "integrator.xml"
state_xml_filename = "state.xml"
state_pdb_filename = "equilibrated.pdb"
system_xml_filename = "system.xml"
checkpoint_filename = "equilibrated.chk"
traj_output_filename = "equilibrated.dcd"
state_data_filename = "state_data.csv"
minimized_pdb_filename = "minimized.pdb"

# Load the AMBER files
print("Creating OpenMM system from AMBER input files...")
prmtop = AmberPrmtopFile(amber_prmtop_file)
inpcrd = AmberInpcrdFile(amber_inpcrd_file)

system = prmtop.createSystem(
    nonbondedMethod=nonbonded_method,
    constraints=constraints,
    temperature=temperature,
    removeCMMotion=remove_cm_motion,
    hydrogenMass=hydrogen_mass,
)

# Fix naked charges in old and new systems
force_dict = {i.__class__.__name__: i for i in system.getForces()}
if 'NonbondedForce' in [i for i in force_dict.keys()]:
    nb_force = force_dict['NonbondedForce']
    for i in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        if sigma == 0*unit.nanometer:
            sigma = 0.6*unit.nanometer
            nb_force.setParticleParameters(i, charge, sigma, epsilon)
        if epsilon == 0*unit.kilojoule_per_mole:
            epsilon = 0.01*unit.kilojoule_per_mole
            nb_force.setParticleParameters(i, charge, sigma, epsilon)

# Add virtual bond between chain R of RBD (CA of VAL172 ) and ACE2 (CA of GLN518)
force = CustomBondForce('0')
force.addBond(2636,8184, [])
system.addForce(force)

# Add a barostat to the system
system.addForce(MonteCarloBarostat(pressure, temperature))

# Make and serialize integrator - Langevin dynamics
print("Serializing integrator to %s" % integrator_xml_filename)
integrator = openmmtools.integrators.LangevinIntegrator(
    temperature,
    collision_rate,  # Friction coefficient
    timestep, 
    constraint_tolerance=1e-5
)
with open(output_prefix + integrator_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(integrator)
    outfile.write(xml)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName("OpenCL")
prop = dict(OpenCLPrecision="mixed")  # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(prmtop.topology, system, integrator, platform, prop)

# Set the particle positions
sim.context.setPositions(inpcrd.positions)

# Minimize the energy
print("Minimising energy...")
print(
    "  initial : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)
sim.minimizeEnergy()
print(
    "  final : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)

# Save minimized structure as a PDB
print("Saving minimized structure  as %s" % minimized_pdb_filename)
with open(output_prefix + minimized_pdb_filename, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True
    )

# set starting velocities:
print("Generating random starting velocities")
sim.context.setVelocitiesToTemperature(temperature)

# write limited state information to standard out:
sim.reporters.append(
    StateDataReporter(
        output_prefix + state_data_filename,
        reportInterval=report_freq,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=nsteps,
        separator="\t",
    )
)

# Write to checkpoint files regularly:
sim.reporters.append(CheckpointReporter(
    file=output_prefix + checkpoint_filename,
    reportInterval=chk_freq
    )
)

# Write out the trajectory
sim.reporters.append(md.reporters.DCDReporter(
    file=output_prefix + traj_output_filename,
    reportInterval=traj_freq
    )
)

# Run NPT dynamics
print("Running dynamics in the NPT ensemble...")
initial_time = time.time()
sim.step(nsteps)
elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = nsteps * timestep
print('    Equilibration took %.3f s for %.3f ns (%8.3f ns/day)' % (elapsed_time / unit.seconds, simulation_time / unit.nanoseconds, simulation_time / elapsed_time * unit.day / unit.nanoseconds))

# Save and serialize the final state
print("Serializing state to %s" % state_xml_filename)
state = sim.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=True
)
with open(output_prefix + state_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(state)
    outfile.write(xml)

# Save the final state as a PDB
print("Saving final state as %s" % state_pdb_filename)
with open(output_prefix + state_pdb_filename, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True
    )

# Save and serialize system
print("Serializing system to %s" % system_xml_filename)
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with open(output_prefix + system_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(system)
    outfile.write(xml)
