#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 30:00 
#
# Set output file
#BSUB -o  repex.soft.10.10.thr.out
#
# Set error file
#BSUB -eo repex.soft.10.10.thr.stderr
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
#BSUB -J "repex 10 10 thr"

source ~/.bashrc
conda activate perses-sims
python run_repex_softened_endstates_flattened_solvent.py "/data/chodera/zhangi/perses_benchmark/neq/10/10/" "vacuum" "THR" 1 5 "forward"

