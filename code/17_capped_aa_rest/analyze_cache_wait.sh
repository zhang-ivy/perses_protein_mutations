#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00 
#
# Set output file
#BSUB -o  11.19.thr.vac.analysis.out
#
# Set error file
#BSUB -eo 11.19.thr.vac.analysis.stderr
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
#BSUB -J "analysis 11 19 thr vac"
#
# wait for another job to be completed
#BSUB -w "done(18988983)"

source ~/.bashrc
conda activate perses-sims
python analyze_cache_backwards.py "/data/chodera/zhangi/perses_benchmark/neq/11/19/" "THR" 1 "vacuum" 5

