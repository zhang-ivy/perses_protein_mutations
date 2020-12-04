#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00
#
# Set output file
#BSUB -o  blues.thr.vacuum.restraints.out
#
# Set error file
#BSUB -eo blues.thr.vacuum.restraint.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "thr vacuum restr"

source ~/.bashrc
conda activate perses-sims
module load cuda/10.1
python generate_blues_cache_aa_restraints.py thr_vacuum.yaml 1 "/data/chodera/zhangi/perses_benchmark/neq/7/21/"
