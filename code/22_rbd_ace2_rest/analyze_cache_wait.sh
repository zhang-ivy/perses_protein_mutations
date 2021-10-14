#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 12:00 
#
# Set output file
#BSUB -o  14.1.analysis.%I.out
#
# Set error file
#BSUB -eo 14.1.analysis.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu ls-gpu lt-gpu lu-gpu lv-gpu lw-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "14.1[1-4]"
# wait for another job to be completed
#BSUB -w "done(79089)"

source ~/.bashrc
#module load cuda/10.1
conda activate perses-rbd-ace2
python analyze_cache_general.py "/data/chodera/zhangi/perses_benchmark/neq/14/1/" "439" "ASN" "LYS" "${LSB_JOBINDEX}"

