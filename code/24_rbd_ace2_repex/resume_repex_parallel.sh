#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 7:00 
#
# Set output file
#BSUB -o  apo_resume.out
#
# Set error file
#BSUB -eo apo_resume.stderr
#
# Specify node group
#BSUB -m "lt-gpu lp-gpu lg-gpu lv-gpu ld-gpu lu-gpu lx-gpu ly-gpu lw-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=2]"
#BSUB -n 4 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "3.apo"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2
export HDF5_USE_FILE_LOCKING=FALSE
build_mpirun_configfile --configfilepath configfile_apo --hostfilepath hostfile_apo "python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/24_rbd_ace2_repex/resume_h_repex.py /data/chodera/zhangi/perses_benchmark/repex/31/3/0/ apo --total_iterations 10000"
mpiexec.hydra -f hostfile_apo -configfile configfile_apo
