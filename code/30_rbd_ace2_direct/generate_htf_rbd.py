import pickle
import os
from perses.app.relative_point_mutation_setup import PointMutationExecutorRBD
from perses.utils.smallmolecules import  render_protein_residue_atom_mapping
from pkg_resources import resource_filename
from simtk import unit
import argparse

# Read args
parser = argparse.ArgumentParser(description='generate htf')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('residue', type=str, help='residue position')
parser.add_argument('mutant', type=str, help='three letter code for amino acid to mutate to')
parser.add_argument('--rbd_file', type=str, help='rbd file')
parser.add_argument('--ace2_file', type=str, help='ace2 file')
args = parser.parse_args()

rbd_file = resource_filename('perses', 'data/rbd-ace2/0_rbd.pdb') if not args.rbd_file else args.rbd_file
ace2_file = resource_filename('perses', 'data/rbd-ace2/0_ace2.pdb') if not args.ace2_file else args.ace2_file

solvent_delivery = PointMutationExecutorRBD(rbd_file,
                        'R',
                        args.residue,
                        args.mutant,
                        ligand_input=ace2_file,
                        ionic_strength=0.15*unit.molar,
                        flatten_torsions=True,
                        flatten_exceptions=True,
                        generate_unmodified_hybrid_topology_factory=True,
                        generate_rest_capable_hybrid_topology_factory=True, 
                        debug_dir=os.path.join(args.dir, "debug/")
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

# Render atom map
atom_map_filename = f'{outdir}/atom_map.png'
render_protein_residue_atom_mapping(apo_htf._topology_proposal, atom_map_filename)

# Save pdbs
traj = md.Trajectory(complex_htf.old_positions(complex_htf.hybrid_positions), md.Topology.from_openmm(complex_htf._topology_proposal.old_topology))
traj.save(os.path.join(outdir, "complex_old.pdb"))

traj = md.Trajectory(complex_htf.new_positions(complex_htf.hybrid_positions), md.Topology.from_openmm(complex_htf._topology_proposal.new_topology))
traj.save(os.path.join(outdir, "complex_new.pdb"))

traj = md.Trajectory(apo_htf.old_positions(apo_htf.hybrid_positions), md.Topology.from_openmm(apo_htf._topology_proposal.old_topology))
traj.save(os.path.join(outdir, "apo_old.pdb"))

traj = md.Trajectory(apo_htf.new_positions(apo_htf.hybrid_positions), md.Topology.from_openmm(apo_htf._topology_proposal.new_topology))
traj.save(os.path.join(outdir, "apo_new.pdb"))

