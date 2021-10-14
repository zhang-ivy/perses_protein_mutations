#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  10.36.ser.analysis.out
#
# Set error file
#BSUB -eo 10.36.ser.analysis.stderr
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
#BSUB -J "analysis 10 36 ser"

source ~/.bashrc
conda activate perses-sims
python analyze_cache.py "/data/chodera/zhangi/perses_benchmark/neq/10/36/" "SER" 0 "solvent" 5

