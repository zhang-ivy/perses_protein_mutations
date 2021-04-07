import os
import yaml

# Load in mutation parameters
stream = open("test.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/24_rbd_ace2_repex/run_repex_parallel.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/repex/31/"

# Define job parameters
wall_time = [24, 12] # hours
memory = [2, 1] # GB
phases = ['complex', 'apo']

for k, v in dictionary.items():
	print(f"prepping dir {k}")
	new = k

	for i, phase in enumerate(phases):

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
			elif "#BSUB -R" in line:
				line = line[:21] + str(memory[i]) + line[23:] 
			elif "#BSUB -J" in line:
				line = line[:10] + str(new) + f'.{phase}"\n'
			elif "mpiexec.hydra" in line:
				hostfile = f"hostfile_{phase}"
				configfile = f"configfile_{phase}"
				line = f"mpiexec.hydra -f {hostfile} -configfile {configfile}"
			elif "build_mpirun_configfile" in line:
				line = f'build_mpirun_configfile --configfilepath configfile_{phase} --hostfilepath hostfile_{phase} "python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/24_rbd_ace2_repex/run_h_repex.py /data/chodera/zhangi/perses_benchmark/repex/31/{new}/0/ {phase}"\n'
			lines_new.append(line)

		# Make dir and save new bash file
		os.system(f"mkdir {os.path.join(out_dir, str(new))}")

		with open(os.path.join(out_dir, str(new), "0" ,f"run_repex_parallel_{phase}.sh"), "w") as f:
			f.writelines(lines_new)
