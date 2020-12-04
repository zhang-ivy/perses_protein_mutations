#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:30 
#
# Set output file
#BSUB -o  aa-neq-2_4s.%I.out
#
# Set error file
#BSUB -eo aa-neq-2_4s.%I.stderr
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
#BSUB -J "2_4s[1-10]"

module load cuda/10.1
python run_neq_distrib_no_ster.py "/data/chodera/zhangi/perses_benchmark/neq/4/2_4/" "solvent" "$((${LSB_JOBINDEX}-1))"

