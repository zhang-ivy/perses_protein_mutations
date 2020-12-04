#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00 
#
# Set output file
#BSUB -o  htf.out
#
# Set error file
#BSUB -eo htf.stderr
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "htf"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python generate_all_htfs.py
