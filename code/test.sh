#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 0:10 
#
# Set output file
#BSUB -o  test.out
#
# Set error file
#BSUB -eo test.stderr
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=1]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "test"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_neq.py "/data/chodera/zhangi/perses_benchmark/neq/0/1/" "solvent"


