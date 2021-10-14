#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00 
#
# Set output file
#BSUB -o  11.19.ala.sol.analysis.out
#
# Set error file
#BSUB -eo 11.19.ala.sol.analysis.stderr
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
#BSUB -J "analysis 11 19 ala sol"

source ~/.bashrc
conda activate perses-sims
python analyze_cache_forward.py "/data/chodera/zhangi/perses_benchmark/neq/11/19/" "ALA" 0 "solvent" 5

