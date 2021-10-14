#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  7.44.thr.analysis.out
#
# Set error file
#BSUB -eo 7.44.thr.analysis.stderr
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
#BSUB -J "analysis 7 44 thr"
#
# wait for another job to be completed
#BSUB -w "done(18590494)"

source ~/.bashrc
conda activate perses-sims
python analyze_cache.py "/data/chodera/zhangi/perses_benchmark/neq/7/44/" "THR" 0 "solvent" 5

