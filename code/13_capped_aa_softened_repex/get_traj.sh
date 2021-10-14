#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  10.34.thr.traj.out
#
# Set error file
#BSUB -eo 10.34.thr.traj.stderr
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
#BSUB -J "traj 10 34 thr"

source ~/.bashrc
conda activate perses-sims
python get_traj_for_state.py "/data/chodera/zhangi/perses_benchmark/neq/10/34/" "solvent" "THR" 0 5 0

