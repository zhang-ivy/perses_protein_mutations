import os
import yaml

# Load in mutation parameters
stream = open("test.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/18_bnbs_rest/run_rest2_template.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/12/"

# Define job parameters
wall_time = 12
memory = 5
script_version = 4

for k, v in dictionary.items():
	new = k
	old, old_aa, new_aa, resid = v

# Edit template bash file
with open(template_file, "r") as f:
    lines_out = f.readlines()

lines_new = []
for line in lines_out:
    if "#BSUB -W" in line:
        line = line[:9] + str(wall_time) + line[11:]
    elif "#BSUB -o" in line:
        line = line[:12] + str(new) + ".%I.out\n"
    elif "#BSUB -eo" in line:
        line = line[:13] + str(new) + ".%I.stderr\n"
    elif "#BSUB -n" in line:
        line = line[:26] + str(memory) + line[27:] 
    elif "#BSUB -J" in line:
        line = line[:13] + str(new) + '[1-4]"\n'
    elif "old=" in line:
        line = f"old={old}\n"
    elif "new=" in line:
        line = f"new={new}\n"
    elif "old_aa=" in line:
        line = f"old_aa={old_aa}\n"
    elif "new_aa=" in line:
        line = f"new_aa={new_aa}\n"
    elif "resid=" in line:
        line = f"resid={resid}\n"
    elif "python generate_rest2_cache" in line:
        line = line[:51] + str(script_version) + line[52:]
    lines_new.append(line)

# Make dir and save new bash file
os.system(f"mkdir {os.path.join(out_dir, str(new))}")

with open(os.path.join(out_dir, str(new), "run_rest2_all.sh"), "w") as f:
     f.writelines(lines_new)
