#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 20:00 
#
# Set output file
#BSUB -o  repex-26c.%I.out
#
# Set error file
#BSUB -eo repex-26c.%I.stderr
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
#BSUB -J "26c[2]"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python run_repex_linear_elec_gp.py "/data/chodera/zhangi/perses_benchmark/repex/26/$((${LSB_JOBINDEX}-1))/" "complex"

