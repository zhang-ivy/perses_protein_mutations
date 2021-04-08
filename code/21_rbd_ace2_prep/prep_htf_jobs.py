import os
import yaml

# Load in mutation parameters
stream = open("test.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/21_rbd_ace2_prep/generate_htf.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/repex/31/"

# Define job parameters
for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k
        old, resid, new_aa = v


		# Edit template bash file
		with open(template_file, "r") as f:
			lines_out = f.readlines()

		lines_new = []
		for line in lines_out:
			if "#BSUB -W" in line:
				line = line[:9] + str(wall_time[i]) + line[10:]
			elif "#BSUB -o" in line:
				line = line[:9] + f" {phase}.out\n"
			elif "#BSUB -eo" in line:
				line = line[:9] + f" {phase}.stderr\n"
			elif "#BSUB -n 1 -R" in line:
				line = line[:26] + str(memory[i]) + line[27:] 
			elif "#BSUB -J" in line:
				line = line[:10] + str(new) + f'.{phase}"\n'
			elif "python " in line:
				line = f"python run_h_repex.py /data/chodera/zhangi/perses_benchmark/repex/31/{new}/0/ {phase}"
			lines_new.append(line)

		# Make dir and save new bash file
		os.system(f"mkdir {os.path.join(out_dir, str(new))}")

		with open(os.path.join(out_dir, str(new), "0" ,f"run_repex_{phase}.sh"), "w") as f:
			f.writelines(lines_new)
