from blues.moves import SideChainMove
from blues.moves import MoveEngine
from blues.simulation import SimulationFactory
import json
from blues.settings import Settings
import openeye.oechem as oechem

# Parse a YAML configuration, return as Dict
cfg = Settings('thr_vacuum.yaml').asDict()
structure = cfg['Structure']

sidechain = SideChainMove(structure, [1], write_move=True)

class SideChainMove2(SideChainMove):
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

sidechain2 = SideChainMove2(structure, [1], [10, 11, 12, 13, 14, 15, 16, 17], write_move=True)
sidechain_mover = MoveEngine(sidechain2)

systems = SystemFactory(structure, sidechain2.atom_indices, cfg['system'])

simulations = SimulationFactory(systems, sidechain_mover, cfg['simulation'], cfg['md_reporters'],
                                cfg['ncmc_reporters'])

outdir = "/data/chodera/zhangi/perses_benchmark/neq/7/13/blues/"

for i in range(100):
    print(f"iter: {i}")
    outfile = os.path.join(outdir, f'thr_{i}.nc')

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
        subprocess.call(['rm', outfile])
        continue

