#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  md.38.0.out
#
# Set error file
#BSUB -eo md.38.0.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=12]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "38 0"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_md_hybrid_flattened.py "/data/chodera/zhangi/perses_benchmark/neq/7/38/" "vacuum" "0" "reverse"

