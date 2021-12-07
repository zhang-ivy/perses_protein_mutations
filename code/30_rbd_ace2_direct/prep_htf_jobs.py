import os
import yaml

# Load in mutation parameters
stream = open("pilot_4.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/30_rbd_ace2_direct/generate_htf.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/14/"

# Define job parameters
for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k
	_, new_aa, resid = v

	# Edit template bash file
	with open(template_file, "r") as f:
		lines_out = f.readlines()

	lines_new = []
	for line in lines_out:
		if "#BSUB -o" in line:
			line = line[:12] + f"{new}.out\n"
		elif "#BSUB -eo" in line:
			line = line[:13] + f"{new}.stderr\n"
		elif "#BSUB -J" in line:
			line = line[:13] + str(new) + '"\n'
		elif "new=" in line:
			line = f"new={new}\n"
		elif "resid=" in line:
			line = f"resid={resid}\n"
		elif "new_aa=" in line:
			line = f"new_aa={new_aa}\n"

		lines_new.append(line)

	# Save bash file
	with open(os.path.join(out_dir, str(new), f"generate_htf.sh"), "w") as f:
		f.writelines(lines_new)
