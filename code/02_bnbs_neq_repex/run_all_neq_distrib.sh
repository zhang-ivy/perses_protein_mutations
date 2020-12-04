#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 20:00 
#
# Set output file
#BSUB -o  8-3c.%I.out
#
# Set error file
#BSUB -eo 8-3c.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=1]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "8-3c[1-100]"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_neq_distrib.py "/data/chodera/zhangi/perses_benchmark/neq/8/3/" "complex" "$((${LSB_JOBINDEX}-1))"

