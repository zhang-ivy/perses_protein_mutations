#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 5:00
#
# Set output file
#BSUB -o  t_repex.18_thr.out
#
# Set error file
#BSUB -eo t_repex.18_thr.stderr
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
#BSUB -J "t_repex_18_thr"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python restart_t_repex_cache_2fs_restraints.py "/data/chodera/zhangi/perses_benchmark/neq/7/18/" 5
