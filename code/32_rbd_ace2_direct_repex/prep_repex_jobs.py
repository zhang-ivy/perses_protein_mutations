import os

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/32_rbd_ace2_direct_repex/run_repex_parallel_complex_scatter.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/repex/32/0/"

# Define job parameters
wall_time = [50] # hours
memory = [7] # GB
phases = ['complex']
n_chunks = 10
name = "32.0"

for i in range(n_chunks):
    print(f"prepping dir {i}")

    lambda_start = i * 0.1
    lambda_end = lambda_start + 0.1

    for i, phase in enumerate(phases):

        # Edit template bash file
        with open(template_file, "r") as f:
            lines_out = f.readlines()

        lines_new = []
        for line in lines_out:
            if "#BSUB -W" in line:
                line = line[:9] + str(wall_time[i]) + line[11:]
            elif "#BSUB -o" in line:
                line = line[:13] + "{phase}.out\n"
            elif "#BSUB -eo" in line:
                line = line[:13] + "{phase}.stderr\n"
            elif "#BSUB -n" in line:
                line = line[:26] + str(memory[i]) + line[27:] 
            elif "#BSUB -J" in line:
                line = line[:13] + name + f'.{i}"\n'
                    elif "new=" in line:
                            line = f"new={i}\n"
            elif "phase=" in line:
                line = f"phase={phase}\n"
            elif "lambda_start=" in line:
                            line = f"lambda_start={lambda_start}\n"
                    elif "lambda_end=" in line:
                            line = f"lambda_end={lambda_end}\n"
            lines_new.append(line)


        with open(os.path.join(out_dir, i, f"run_repex_parallel_{phase}_scatter.sh"), "w") as f:
            f.writelines(lines_new)
