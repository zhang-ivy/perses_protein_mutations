#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Set output file
#BSUB -o  T42A.2.5a.%I.out
#
# Set error file
#BSUB -eo T42A.2.5a.%I.stderr
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
#BSUB -J "T42A25a[22-24]"

module load cuda/10.1
python run_repex_shorter.py "/data/chodera/zhangi/perses_benchmark/repex/"${LSB_JOBINDEX}"/0/" "apo"

