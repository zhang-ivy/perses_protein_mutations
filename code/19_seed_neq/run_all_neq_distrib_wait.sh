#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o  13.73.apo.%I.out
#
# Set error file
#BSUB -eo 13.73.apo.%I.stderr
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
#BSUB -J "13.73.apo[1-100]"
# wait for another job to be completed
#BSUB -w "done(19948233)"

source ~/.bashrc
conda activate perses-sims-oct2020
module load cuda/10.1
python run_neq_distrib_seed_general.py "/data/chodera/zhangi/perses_benchmark/neq/13/73/" "apo" "$((${LSB_JOBINDEX}-1))" "ala" "asp" 1 # lowercase
