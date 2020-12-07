#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  10.11.solvent.%I.out
#
# Set error file
#BSUB -eo 10.11.solvent.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=3]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "10.11.solvent[1-100]"

source ~/.bashrc
conda activate perses-sims-oct2020
module load cuda/10.1
python run_neq_distrib_longer_flattened.py "/data/chodera/zhangi/perses_benchmark/neq/10/11/" "solvent" "$((${LSB_JOBINDEX}-1))" 1 1 
