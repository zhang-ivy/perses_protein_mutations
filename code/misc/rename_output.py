import os
import subprocess

dir_0 = "/data/chodera/zhangi/perses_benchmark/neq/6/0/"
dir_1 = "/data/chodera/zhangi/perses_benchmark/neq/6/1/"
# for file in os.listdir(dir_0):
#     if file.startswith("0_apo"):
file_ends = ['forward_eq_old.pdb', 'forward_neq_new.pdb', 'forward_neq_old.pdb', 'forward.npy', 
'reverse_eq_new.pdb', 'reverse_neq_new.pdb', 'reverse_neq_old.pdb', 'reverse.npy']

# counter = 70
# for i in range(10):
#     for j in range(3):
#         for file_end in file_ends:
#             old_file = os.path.join(dir_1, f"1_apo_{i}_{j}_{file_end}")
#             new_file = f"1_apo_{counter}_" + file_end
#             print(old_file, " ", new_file)
#             command = ["mv", old_file, os.path.join(dir_1, new_file)]
#             try:
#                 subprocess.run(command, check=True)
#             except Exception as e:
#                 print(f"Failed: {e.stdout}") # Note the stdout is empty
#         counter += 1

counter = 90
for i in range(10):
    j = 0
    for file_end in file_ends:
        old_file = os.path.join(dir_0, f"0_complex_{i}_{j}_{file_end}")
        new_file = f"0_complex_{counter}_" + file_end
        print(old_file, " ", new_file)
        command = ["mv", old_file, os.path.join(dir_0, new_file)]
        try:
            subprocess.run(command, check=True)
        except Exception as e:
            print(f"Failed: {e.stdout}") # Note the stdout is empty
    counter += 1