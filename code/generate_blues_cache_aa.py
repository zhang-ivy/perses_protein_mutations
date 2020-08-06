from blues.moves import SideChainMove
from blues.moves import MoveEngine
from blues.simulation import SimulationFactory, SystemFactory, BLUESSimulation
import json
from blues.settings import Settings
import openeye.oechem as oechem
import os
import subprocess
from blues.reporters import NetCDF4Reporter
import argparse
from simtk import openmm

import pickle
import numpy as np
import os

# Read args
parser = argparse.ArgumentParser(description='run perses protein mutation on capped amino acid')
parser.add_argument('yaml', type=str, help='yaml file for blues')
parser.add_argument('resid', type=int, help='resid of sidechain that needs to be sampled')
parser.add_argument('outdir', type=str, help='outdir for blues nc files')
args = parser.parse_args()

# Read in htf
i = os.path.basename(os.path.dirname(args.outdir))
with open(os.path.join(args.outdir, f"{i}_vacuum.pickle"), 'rb') as f:
    htf = pickle.load(f)

# Read in PDB
from simtk.openmm import app
pdb = app.PDBFile(os.path.join(outdir, "blues/thr_vacuum.pdb"))

# Set periodic box vectors
vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()
pdb.topology.setPeriodicBoxVectors(vectors)

# Parse a YAML configuration, return as Dict
cfg = Settings(os.path.join(outdir, "blues/thr_vacuum.yaml")).asDict()

# Overwrite the parmed structure object in the cfg dict with pdb object
cfg['Structure'] = pdb

# Initialize structure variable to be used in later BLUES calls
structure = cfg['Structure']

# Create OpenMM system generator
from openmmforcefields.generators import SystemGenerator
import simtk.unit as unit
from simtk.openmm import app
temperature = 300 * unit.kelvin
forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
# forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 0.00025, 'constraints' : app.HBonds, 'hydrogenMass' : 4 * unit.amus}
forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 0.00025, 'hydrogenMass' : 4 * unit.amus}
periodic_forcefield_kwargs = None
nonperiodic_forcefield_kwargs = {'nonbondedMethod': app.NoCutoff}
small_molecule_forcefields = 'gaff-2.11'
molecules = []
system_generator = SystemGenerator(forcefields=forcefield_files,
                                        forcefield_kwargs=forcefield_kwargs,
                                        periodic_forcefield_kwargs=periodic_forcefield_kwargs,
                                        nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs,
                                        small_molecule_forcefield=small_molecule_forcefields,
                                        molecules=molecules,
                                        cache=None)

## Subclass the BLUES SideChainMove object to avoid using parmed.Structure object 
## and handle an OpenMM modeller object instead (Changes involve how model positions are set)

class SideChainMoveOpenMM2(SideChainMove):  
    def __init__(self, structure, residue_list, qry_indices, verbose=False, write_move=False):
            self.structure = structure
            self.molecule = self._pmdStructureToOEMol()
            self.residue_list = residue_list
            self.all_atoms = [atom.index for atom in self.structure.topology.atoms()]
            self.rot_atoms, self.rot_bonds, self.qry_atoms = self.getRotBondAtoms(qry_indices)
            self.atom_indices = self.rot_atoms
            self.verbose = verbose
            self.write_move = write_move

    def getTargetAtoms(self, molecule, residue_list, qry_indices,):
            """Takes an OpenEye molecule and a list of residue numbers then
            generates a dictionary containing all the atom pointers and indicies for the
            non-backbone, atoms of those target residues, as well as a list of backbone atoms.
            Note: The atom indicies start at 0 and are thus -1 from the PDB file indicies
            Parameters
            ----------
            molecule : oechem.OEMolecule
                The OEmolecule of the simulated system.
            qry_indices : list of int
                List containing the atom indices of the target atoms.
            residue_list : list of int
                List containing the residue numbers of the sidechains to be rotated.
            Returns
            -------
            backbone_atoms : list of int
                List containing the atom indices of the backbone atoms to be rotated.
            qry_atoms : dict of oechem.OEAtomBase
                Dictionary containing all the atom pointers (as OpenEye objects) that
                make up the given residues.
            """

            # create and clear dictionary to store atoms that make up residue list
            qry_atoms = {}
            qry_atoms.clear()

            reslib = []

            #print('Searching residue list for atoms...')
            # loop through all the atoms in the PDB OEGraphMol structure
            for atom in molecule.GetAtoms():
                # check if the atom is in backbone
                if atom.GetIdx() in qry_indices:
                    # if heavy, find what residue it is associated with
                    myres = oechem.OEAtomGetResidue(atom)
                    # check if the residue number is amongst the list of residues
                    if myres.GetResidueNumber() in residue_list and myres.GetName() != "HOH":
                        # store the atom location in a query atom dict keyed by its atom index
                        qry_atoms.update({atom: atom.GetIdx()})
                        #print('Found atom %s in residue number %i %s'%(atom,myres.GetResidueNumber(),myres.GetName()))
                        if myres not in reslib:
                            reslib.append(myres)
            return qry_atoms
        
    def getRotBondAtoms(self, qry_indices):
        """This function is called on class initialization.
        Takes in a PDB filename (as a string) and list of residue numbers.  Returns
        a nested dictionary of rotatable bonds (containing only heavy atoms), that are keyed by residue number,
        then keyed by bond pointer, containing values of atom indicies [axis1, axis2, atoms to be rotated]
        Note: The atom indicies start at 0, and are offset by -1 from the PDB file indicies
        Returns
        -------
        rot_atoms : dict
            Dictionary of residues, bonds and atoms to be rotated
        rot_bonds : dict of oechem.OEBondBase
            Dictionary containing the bond pointers of the rotatable bonds.
        qry_atoms : dict of oechem.OEAtomBase
            Dictionary containing all the atom pointers (as OpenEye objects) that
            make up the given residues.
        """
        backbone_atoms = [idx for idx in self.all_atoms if idx not in qry_indices]
        
        # Generate dictionary containing locations and indicies of heavy residue atoms
        #print('Dictionary of all query atoms generated from residue list\n')
        qry_atoms = self.getTargetAtoms(self.molecule, self.residue_list, qry_indices)

        # Identify bonds containing query atoms and return dictionary of indicies
        rot_bonds = self.findHeavyRotBonds(self.molecule, qry_atoms)

        # Generate dictionary of residues, bonds and atoms to be rotated
        rot_atoms = self.getRotAtoms(rot_bonds, self.molecule, backbone_atoms)
        return rot_atoms, rot_bonds, qry_atoms

    def move(self, context, verbose=False):
        """Rotates the target atoms around a selected bond by angle theta and updates
        the atom coordinates in the parmed structure as well as the ncmc context object
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        verbose : bool, default=False
            Enable verbosity to print out detailed information of the rotation.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose positions were changed by this function.
        """

        import copy
        import numpy
        from simtk import unit
    
        # Determine the axis, theta, residue, and bond + atoms to be rotated
        theta, target_atoms, res, bond = self.chooseBondandTheta()
        print('Rotating bond: %s in resnum: %s by %.2f radians' % (bond, res, theta))

        # Retrieve the current positions
        initial_positions = context.getState(getPositions=True).getPositions(asNumpy=True)
        nc_positions = copy.deepcopy(initial_positions)

        model = copy.copy(self.structure)
        atoms = list(model.topology.atoms())

        # Set the modeller model to the same coordinates as the context
        for idx, atom in enumerate(self.all_atoms):
            if self.verbose:
                print('Before:')
                print(atom, idx)
                print(nc_positions[atom], model.positions[atom])

            model.positions[atom] = nc_positions[atom]

            if self.verbose:
                print('After:')
                print(nc_positions[atom], model.positions[atom])

        positions = model.positions

        # Find the rotation axis using the updated positions
        axis1 = target_atoms[0]
        axis2 = target_atoms[1]
        rot_axis = (positions[axis1] - positions[axis2]) / positions.unit

        # Calculate the rotation matrix
        rot_matrix = self.rotation_matrix(rot_axis, theta)

        # Apply the rotation matrix to the target atoms
        for idx, atom in enumerate(target_atoms):

            my_position = positions[atom]

            if self.verbose:
                print('The current position for %i is: %s' % (atom, my_position))

            # Find the reduced position (substract out axis)
            red_position = (my_position - model.positions[axis2])._value
            # Find the new positions by multiplying by rot matrix
            new_position = numpy.dot(rot_matrix, red_position) * positions.unit + positions[axis2]

            if self.verbose: print("The new position should be:", new_position)

            positions[atom] = new_position
            # Update the parmed model with the new positions
            model.positions[atom][0] = new_position[0] 
            model.positions[atom][1] = new_position[1] 
            model.positions[atom][2] = new_position[2]

            # Update the copied ncmc context array with the new positions
            nc_positions[atom][0] = model.positions[atom][0] 
            nc_positions[atom][1] = model.positions[atom][1] 
            nc_positions[atom][2] = model.positions[atom][2]
  
            if self.verbose:
                print('The updated position for this atom is:', model.positions[atom])

        # update the actual ncmc context object with the new positions
        context.setPositions(nc_positions)

        # update the class structure positions
        self.structure.positions = model.positions

#         if self.write_move:
#             filename = 'sc_move_%s_%s_%s.pdb' % (res, axis1, axis2)
#             mod_prot = model.save(filename, overwrite=True)
        return context

sidechain = SideChainMoveOpenMM2(structure, [1], [10, 11, 12, 13, 14, 15, 16, 17], write_move=True)
sidechain_mover = MoveEngine(sidechain)

## Subclass the BLUES SystemFactory object to take in an OpenMM SystemGenerator

class SystemFactoryOpenMM(SystemFactory):
    def __init__(self, structure, atom_indices, system_generator, config=None):
        self.structure = structure
        self.atom_indices = atom_indices
        self._config = config
    
        # If parameters for generating the openmm.System is given, make them.
        if self._config:
            if 'alchemical' in self._config.keys():
                self.alch_config = self._config.pop('alchemical')
            else:
                # Use function defaults if none is provided
                self.alch_config = {}
            self.md =  SystemFactoryOpenMM.generateSystem(self, self.structure, system_generator, **self._config)
            self.alch = SystemFactory.generateAlchSystem(self.md, self.atom_indices, **self.alch_config)
              
    def generateSystem(self, structure, system_generator, **kwargs):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file. This function is just a wrapper for parmed Structure.createSystem().
        Parameters
        ----------
        structure : app.PDBFile
            The app.PDBFile of the molecular system to be simulated
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            simtk.openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        constraints : None, app.HBonds, app.HAngles, or app.AllBonds
            Which type of constraints to add to the system (e.g., SHAKE). None
            means no bonds are constrained. HBonds means bonds with hydrogen are
            constrained
        rigidWater : bool=True
            If True, water is kept rigid regardless of the value of constraints.
            A value of False is ignored if constraints is not None.
        implicitSolvent : None, app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2
            The Generalized Born implicit solvent model to use.
        implicitSolventKappa : float or 1/distance Quantity = None
            This is the Debye kappa property related to modeling saltwater
            conditions in GB. It should have units of 1/distance (1/nanometers
            is assumed if no units present). A value of None means that kappa
            will be calculated from implicitSolventSaltConc (below)
        implicitSolventSaltConc : float or amount/volume Quantity=0 moles/liter
            If implicitSolventKappa is None, the kappa will be computed from the
            salt concentration. It should have units compatible with mol/L
        temperature : float or temperature Quantity = 298.15 kelvin
            This is only used to compute kappa from implicitSolventSaltConc
        soluteDielectric : float=1.0
            The dielectric constant of the protein interior used in GB
        solventDielectric : float=78.5
            The dielectric constant of the water used in GB
        useSASA : bool=False
            If True, use the ACE non-polar solvation model. Otherwise, use no
            SASA-based nonpolar solvation model.
        removeCMMotion : bool=True
            If True, the center-of-mass motion will be removed periodically
            during the simulation. If False, it will not.
        hydrogenMass : float or mass quantity = None
            If not None, hydrogen masses will be changed to this mass and the
            difference subtracted from the attached heavy atom (hydrogen mass
            repartitioning)
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        flexibleConstraints : bool=True
            If False, the energies and forces from the constrained degrees of
            freedom will NOT be computed. If True, they will (but those degrees
            of freedom will *still* be constrained).
        verbose : bool=False
            If True, the progress of this subroutine will be printed to stdout
        splitDihedrals : bool=False
            If True, the dihedrals will be split into two forces -- proper and
            impropers. This is primarily useful for debugging torsion parameter
            assignments.
        Returns
        -------
        openmm.System
            System formatted according to the PDB file
        Notes
        -----
        This function calls prune_empty_terms if any Topology lists have
        changed.
        """
        
        return system_generator.create_system(structure.topology)

# Instantiate (modified) BLUES SystemFactory
systems = SystemFactoryOpenMM(structure, sidechain.atom_indices, system_generator, cfg['system'])

## Subclass BLUES SimulationFactory object to avoid using parmed.Structure object, and instead use an
## OpenMM Modeller object. (Changes involves how box vectors are checked/set)

class SimulationFactoryOpenMM(SimulationFactory):
    
    @classmethod
    def generateSimFromStruct(cls, structure, system, integrator, platform=None, properties={}, **kwargs):
        """Generate the OpenMM Simulation objects from a given parmed.Structure()
        Parameters
        ----------
        structure : parmed.Structure
            ParmEd Structure object of the entire system to be simulated.
        system : openmm.System
            The OpenMM System object corresponding to the reference system.
        integrator : openmm.Integrator
            The OpenMM Integrator object for the simulation.
        platform : str, default = None
            Valid choices: 'Auto', 'OpenCL', 'CUDA'
            If None is specified, the fastest available platform will be used.
        Returns
        -------
        simulation : openmm.Simulation
            The generated OpenMM Simulation from the parmed.Structure, openmm.System,
            amd the integrator.
        """
        #Specifying platform properties here used for local development.
        if platform is None:
            #Use the fastest available platform
            simulation = app.Simulation(structure.topology, system, integrator)
        else:
            platform = openmm.Platform.getPlatformByName(platform)
            #Make sure key/values are strings
            properties = {str(k): str(v) for k, v in properties.items()}
            simulation = app.Simulation(structure.topology, system, integrator, platform, properties)

        # Set initial positions/velocities
        if structure.topology.getPeriodicBoxVectors():
            simulation.context.setPeriodicBoxVectors(*structure.topology.getPeriodicBoxVectors())
        simulation.context.setPositions(structure.positions)
        simulation.context.setVelocitiesToTemperature(integrator.getTemperature())

        return simulation

# simulations = SimulationFactory(systems, sidechain_mover, cfg['simulation'], cfg['md_reporters'],
#                                 cfg['ncmc_reporters'])

# Instantiate BLUES SimulationFactory
simulations = SimulationFactoryOpenMM(systems, sidechain_mover, cfg['simulation'], cfg['md_reporters'],
                                cfg['ncmc_reporters'])

for i in range(100):
    print(f"iter: {i}")
    outfile = os.path.join(args.outdir, f'blues/thr_{i}.nc')

    # Manually set new reporter for each iteration i
    reporter = NetCDF4Reporter(outfile, reportInterval=500)
    cfg['md_reporters'] = [cfg['md_reporters'][0]] + [reporter] + cfg['md_reporters'][2:] 

    # Run simulation
    simulations = SimulationFactory(systems, sidechain_mover, cfg['simulation'], cfg['md_reporters'],
                                cfg['ncmc_reporters'])
    simulations.md.minimizeEnergy(maxIterations=0)
    simulations.md.step(500)
    blues = BLUESSimulation(simulations, cfg['simulation'])
    try:
        blues.run()
    except:
        # subprocess.call(['rm', outfile])
        continue

