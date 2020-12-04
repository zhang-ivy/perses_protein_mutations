#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00
#
# Set output file
#BSUB -o  aa-eq.6-no-alc.out
#
# Set error file
#BSUB -eo aa-eq.6-no-alc.stderr
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=15]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "6-eq-no-alc"

python run_equilibration.py "/data/chodera/zhangi/perses_benchmark/neq/3/6/" "solvent"
