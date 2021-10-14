#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 12:00
#
# Set output file
#BSUB -o %I.out
#
# Set error file
#BSUB -eo %I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu lt-gpu lu-gpu lv-gpu lw-gpu lx-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "12.99[1-4]"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
