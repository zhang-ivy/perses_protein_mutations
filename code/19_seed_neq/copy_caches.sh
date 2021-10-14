#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 0:30 
#
# Set output file
#BSUB -o  13.101.copy.out
#
# Set error file
#BSUB -eo 13.101.copy.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=2]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "13.101.c"
# wait for another job to be completed
#BSUB -w "done(79102)"

source ~/.bashrc
cp /data/chodera/zhangi/perses_benchmark/neq/12/101/*.npy /data/chodera/zhangi/perses_benchmark/neq/13/101/
