import os
import yaml

# Load in mutation parameters
stream = open("pilot_1.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/29_bnbs_reprepped/run_neq_vanilla.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/23/"

# Define job parameters
wall_time = [5, 3] # hours
memory = [7, 3] # GB
eq_length = 2 # ns
neq_length = 1 # ns
phases = ['complex', 'apo']

for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k
	old_aa, new_aa, resid, _ = v

	for i, phase in enumerate(phases):

		# Edit template bash file
		with open(template_file, "r") as f:
			lines_out = f.readlines()

		lines_new = []
		for line in lines_out:
			if "#BSUB -W" in line:
				line = line[:9] + str(wall_time[i]) + line[10:]
			elif "#BSUB -o" in line:
				line = line[:13] + str(new) + f".{phase}.%I.out\n"
			elif "#BSUB -eo" in line:
				line = line[:13] + str(new) + f".{phase}.%I.stderr\n"
			elif "#BSUB -n" in line:
				line = line[:26] + str(memory[i]) + line[27:] 
			elif "#BSUB -J" in line:
				line = line[:13] + str(new) + f'.{phase}.[1-100]"\n'
			elif "#BSUB -w" in line:
				job = job_apo if phase == 'apo' else job_complex
				line = line[:9] + f'"done({job})"\n'
			elif "new=" in line:
				line = f"new={new}\n"
			elif "old_aa=" in line:
				line = f"old_aa={old_aa.lower()}\n"
			elif "new_aa=" in line:
				line = f"new_aa={new_aa.lower()}\n"
			elif "phase=" in line:
				line = f"phase={phase}\n"
			elif line.startswith("eq_length="):
				line = f"eq_length={eq_length}\n"
			elif line.startswith("neq_lengt="):
				line = f"neq_length={neq_length}\n"
			elif "resid=" in line:
				line = f"resid={resid}\n"
			#elif "cp" in line:
			#	if os.path.exists(os.path.join(out_dir, f"{new}_{phase}.pickle")):
			#		continue
			lines_new.append(line)


		with open(os.path.join(out_dir, str(new), f"run_neq_vanilla_{phase}.sh"), "w") as f:
			f.writelines(lines_new)
