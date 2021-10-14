#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  12.47.analysis.%I.out
#
# Set error file
#BSUB -eo 12.47.analysis.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu ls-gpu lt-gpu lu-gpu lv-gpu lw-gpu lx-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "12.47*[1-4]"
#
# wait for another job to be completed
#BSUB -w "done(19217124)"

source ~/.bashrc
conda activate perses-sims-oct2020
python analyze_cache_general.py "/data/chodera/zhangi/perses_benchmark/neq/12/47/" "29" "PHE" "TYR" "${LSB_JOBINDEX}"
