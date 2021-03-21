import os
import yaml

# Load in mutation parameters
stream = open("test.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/22_rbd_ace2_rest/run_rest2_template.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/14/"

# Define job parameters
phases =  ['complex', 'apo']
#wall_time = [30, 15]
#memory = [10, 5]
wall_time = [2, 2]
memory = [10, 5]
script_version = 4

for k, v in dictionary.items():
	new = k
	old, old_aa, new_aa, resid, job = v
	
	for i, phase in enumerate(phases):

		# Edit template bash file
		with open(template_file, "r") as f:
		    lines_out = f.readlines()

		lines_new = []
		for line in lines_out:
		    if "#BSUB -W" in line:
		        line = line[:9] + str(wall_time[i]) + line[11:]
		    elif "#BSUB -o" in line:
		        line = line[:12] + str(new) + ".%I.out\n"
		    elif "#BSUB -eo" in line:
		        line = line[:13] + str(new) + ".%I.stderr\n"
		    elif "#BSUB -n" in line:
		        line = line[:26] + str(memory[i]) + line[28:] 
		    elif "#BSUB -J" in line:
		        job_range = "1-2" if phase == 'apo' else "3-4"
		        line = line[:13] + str(new) + f'[{job_range}]"\n'
		    elif "#BSUB -w" in line:
	            line = line[:9] + f'"done({job})"\n'
		    elif "new=" in line:
		        line = f"new={new}\n"
		    elif "old=" in line:
			    line = f"old={old}\n"
		    elif "oldaa=" in line:
		        line = f"old_aa={old_aa}\n"
		    elif "new_aa=" in line:
		        line = f"new_aa={new_aa}\n"
		    elif "resid=" in line:
		        line = f"resid={resid}\n"
		    lines_new.append(line)

		# Make dir and save new bash file
		os.system(f"mkdir {os.path.join(out_dir, str(new))}")

		with open(os.path.join(out_dir, str(new), f"run_rest2_{phase}.sh"), "w") as f:
		     f.writelines(lines_new)
