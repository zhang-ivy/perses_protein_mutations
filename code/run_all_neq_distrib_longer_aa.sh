#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 15:00 
#
# Set output file
#BSUB -o  neq-4.1_12v.%I.out
#
# Set error file
#BSUB -eo neq-4.1_12v.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "1-12v[1-100]"

module load cuda/10.1
python run_neq_distrib_longer_aa.py "/data/chodera/zhangi/perses_benchmark/neq/4/1_12/" "vacuum" "$((${LSB_JOBINDEX}-1))"
