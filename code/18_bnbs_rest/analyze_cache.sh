#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 12:00 
#
# Set output file
#BSUB -o  12.39.ala.complex.analysis.out
#
# Set error file
#BSUB -eo 12.39.ala.complex.analysis.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu ls-gpu lt-gpu lu-gpu lv-gpu lw-gpu lx-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=20]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "analysis 12 39 ala com*"
#

source ~/.bashrc
conda activate perses-sims-oct2020
python analyze_cache_forward.py "/data/chodera/zhangi/perses_benchmark/neq/12/39/" "ALA" 0 "complex" 5

