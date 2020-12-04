#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o 7.29.%I.out
#
# Set error file
#BSUB -eo 7.29.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=2]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "7.29[1-100]"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_neq_distrib_cache_2fs_dominic.py "/data/chodera/zhangi/perses_benchmark/neq/7/29/" "vacuum" "$((${LSB_JOBINDEX}-1))" "ala" "thr"
