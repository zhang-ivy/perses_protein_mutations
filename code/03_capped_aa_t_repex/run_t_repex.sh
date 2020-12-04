#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00
#
# Set output file
#BSUB -o  t_repex.39_ala.out
#
# Set error file
#BSUB -eo t_repex.39_ala.stderr
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
#BSUB -J "t_repex_39_ala"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python generate_t_repex_cache_2fs_flattened.py "/data/chodera/zhangi/perses_benchmark/neq/7/39/" "vacuum" "ALA" 0 5 "forward"
