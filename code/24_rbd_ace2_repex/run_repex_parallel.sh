#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o  complex.out
#
# Set error file
#BSUB -eo complex.stderr
#
# Specify node group
#BSUB -m "lt-gpu lp-gpu lg-gpu lv-gpu ld-gpu lu-gpu lx-gpu ly-gpu lw-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=20]"
#BSUB -n 4 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "3.0.com"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2
build_mpirun_configfile --configfilepath configfile --hostfilepath hostfile "python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/24_rbd_ace2_repex/run_h_repex.py /data/chodera/zhangi/perses_benchmark/repex/31/3/0/ complex"
mpiexec.hydra -f hostfile -configfile configfile

