import os
import yaml

# Load in mutation parameters
stream = open("pilot_5.yaml", 'r')
dictionary = yaml.load(stream)

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/33_atom_mapping/run_neq_template.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/neq/15/"

# Define job parameters
wall_time = [6, 3] # hours
memory = [7, 3] # GB
eq_length = 1 # ns
neq_length = 1 # ns
chain_A = 0
chain_B = 2
phases = ['complex', 'apo']

for k, v in dictionary.items():
    print(f"prepping dir {k}")
    new = k
    old_aa, new_aa, resid = v

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
            elif "phase=" in line:
                line = f"phase={phase}\n"
            elif "eq_length=" in line:
                line = f"eq_length={eq_length}\n"
            elif "neq_length=" in line:
                line = f"neq_length={neq_length}\n"
            elif "chain_A=" in line:
                line = f"chain_A={chain_A}\n"
            elif "chain_B=" in line:
                line = f"chain_B={chain_B}\n"
            lines_new.append(line)

        with open(os.path.join(out_dir, str(new), f"run_neq_{phase}.sh"), "w") as f:
                f.writelines(lines_new)
