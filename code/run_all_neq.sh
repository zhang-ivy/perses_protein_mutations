#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  aa-0s-neq.%I.out
#
# Set error file
#BSUB -eo aa-0s-neq.%I.stderr
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
#BSUB -J "aa-0s-neq[1-12]"

module load cuda/10.1
python run_neq.py "/data/chodera/zhangi/perses_benchmark/neq/0/$((${LSB_JOBINDEX}-1))/" "solvent"


