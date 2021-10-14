#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  16.3.solvent.%I.out
#
# Set error file
#BSUB -eo 16.3.solvent.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lu-gpu lv-gpu ld-gpu lt-gpu lx-gpu ly-gpu lp-gpu lg-gpu lw-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=3]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "16.3.sol[1-100]"

source ~/.bashrc
conda activate perses-counterion
python run_neq_distrib_flattened.py "/data/chodera/zhangi/perses_benchmark/neq/16/3/" "solvent" "$((${LSB_JOBINDEX}-1))" 1 1 
