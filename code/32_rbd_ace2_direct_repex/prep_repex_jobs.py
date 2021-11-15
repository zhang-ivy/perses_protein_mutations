import os

# Define file paths
template_file = "/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/32_rbd_ace2_direct_repex/run_repex_parallel_complex.sh"
out_dir = "/data/chodera/zhangi/perses_benchmark/repex/32/1/"

# Define job parameters
wall_time = [50] # hours
memory = [20] # GB
phases = ['complex']
n_chunks = 10
name = "32.1"

for chunk in range(n_chunks):
    print(f"prepping dir {chunk}")

    lambda_start = chunk * 0.1
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
                line = line[:10] + f"{phase}.out\n"
            elif "#BSUB -eo" in line:
                line = line[:10] + f"{phase}.stderr\n"
            elif "#BSUB -R" in line:
                line = line[:21] + str(memory[i]) + line[23:] 
            elif "#BSUB -J" in line:
                line = line[:10] + name + f'.{chunk}.{phase}"\n'
            elif "new=" in line:
                line = f"new={chunk}\n"
            elif "phase=" in line:
                line = f"phase={phase}\n"
            elif "lambda_start=" in line:
                line = f"lambda_start={lambda_start}\n"
            elif "lambda_end=" in line:
                line = f"lambda_end={lambda_end}\n"
            lines_new.append(line)


        with open(os.path.join(out_dir, str(chunk), f"run_repex_parallel_{phase}.sh"), "w") as f:
            f.writelines(lines_new)
