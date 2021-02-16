#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00 
#
# Set output file
#BSUB -o  12.73.analysis.%I.out
#
# Set error file
#BSUB -eo 12.73.analysis.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu ls-gpu lt-gpu lu-gpu lv-gpu lw-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "12.73*[1-4]"
# wait for another job to be completed
#BSUB -w "done(19948228)"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
python analyze_cache_general.py "/data/chodera/zhangi/perses_benchmark/neq/12/73/" "39" "ALA" "ASP" "${LSB_JOBINDEX}"

