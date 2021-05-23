import os
import yaml

# Load in mutation parameters
stream = open("neq_2.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/27_bnbs_counterion/run_all_neq_distrib_wait_template.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/18/"

# Define job parameters
wall_time = [5, 2] # hours
memory = [7, 3] # GB
length = 1 # ns
phases = ['complex', 'apo']

for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k
	resid, old_aa, new_aa, job_apo, job_complex = v

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
			elif "length=" in line:
				line = f"length={length}\n"
			elif "resid=" in line:
				line = f"resid={resid}\n"
			lines_new.append(line)


		with open(os.path.join(out_dir, str(new), f"run_all_neq_distrib_wait_{phase}.sh"), "w") as f:
			f.writelines(lines_new)
