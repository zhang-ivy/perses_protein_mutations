#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 30:00 
#
# Set output file
#BSUB -o  repex.soft.5.1.out
#
# Set error file
#BSUB -eo repex.soft.5.1.stderr
#
# Specify node group
#BSUB -q cpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#
# job name (default = name of script file)
#BSUB -J "repex soft 5 1"

source ~/.bashrc
conda activate perses-sims
python run_repex_softened_endstates_thr_ala.py "/data/chodera/zhangi/perses_benchmark/neq/10/5/" 1 5

