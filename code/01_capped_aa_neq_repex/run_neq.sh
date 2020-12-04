#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00
#
# Set output file
#BSUB -o  aa-neq-2_1v.out
#
# Set error file
#BSUB -eo aa-neq-2_1v.stderr
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "aa-neq-2_1v"

python run_neq_ala_thr.py "/data/chodera/zhangi/perses_benchmark/neq/4/2_1/" "vacuum"
