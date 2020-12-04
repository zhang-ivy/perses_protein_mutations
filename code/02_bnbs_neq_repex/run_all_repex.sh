#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 4:00 
#
# Set output file
#BSUB -o  repex-29v.%I.out
#
# Set error file
#BSUB -eo repex-29v.%I.stderr
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
#BSUB -J "repex_29v[3]"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_repex_linear_elec.py "/data/chodera/zhangi/perses_benchmark/repex/29/$((${LSB_JOBINDEX}-1))/" "vacuum"

