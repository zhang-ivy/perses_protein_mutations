#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  12.20.complex.%I.out
#
# Set error file
#BSUB -eo 12.20.complex.%I.stderr
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
#BSUB -J "12.20.com[1-100]"

source ~/.bashrc
conda activate perses-sims-oct2020
module load cuda/10.1
python run_neq_distrib_rest2.py "/data/chodera/zhangi/perses_benchmark/neq/12/20/" "complex" "$((${LSB_JOBINDEX}-1))" "ALA" "THR" 1 
