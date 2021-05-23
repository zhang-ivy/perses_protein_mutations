import os
import yaml

# Load in mutation parameters
stream = open("test.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/27_bnbs_counterion/generate_htfs.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/17/"

barnase_path = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/mmc2_barnase.pdb"
barstar_path = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/mmc2_barstar.pdb"

# Define job parameters
for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k
	resid, new_aa, protein = v

	if protein == "barnase":
		protein_path = barnase_path
		ligand_path = barstar_path
	elif protein == "barstar":
		protein_path = barstar_path
		ligand_path = barnase_path
	else:
		raise Exception(f"Invalid protein name in yaml. You specified {protein} but valid options are: 'barnase' or 'barstar'")

	# Edit template bash file
	with open(template_file, "r") as f:
		lines_out = f.readlines()

	lines_new = []
	for line in lines_out:
		if "#BSUB -o" in line:
			line = line[:13] + f"{new}.out\n"
		elif "#BSUB -eo" in line:
			line = line[:13] + f"{new}.stderr\n"
		elif "#BSUB -J" in line:
			line = line[:13] + str(new) + '"\n'
		elif "dir=" in line:
			line = f"dir={new}\n"
		elif "resid=" in line:
			line = f"resid={resid}\n"
		elif "new_aa=" in line:
			line = f"new_aa={new_aa}\n"
		elif "protein_path=" in line:
			line = f"protein_path={protein_path}\n"
		elif "ligand_path=" in line:
			line = f"ligand_path={ligand_path}\n"

		lines_new.append(line)

	# Save bash file
	with open(os.path.join(out_dir, str(new), f"generate_htfs.sh"), "w") as f:
		f.writelines(lines_new)
