#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  11.19.thr.sol.out
#
# Set error file
#BSUB -eo 11.19.thr.sol.stderr
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
#BSUB -J "11 19 thr sol"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
python generate_rest2_cache.py "/data/chodera/zhangi/perses_benchmark/neq/11/19/" "solvent" "THR" 1 5 1200 "backward"

