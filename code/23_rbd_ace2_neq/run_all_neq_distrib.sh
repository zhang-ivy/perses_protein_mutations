#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  15.1.apo.%I.out
#
# Set error file
#BSUB -eo 15.1.apo.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "15.1.apo[1-100]"

source ~/.bashrc
conda activate perses-rbd-ace2
#module load cuda/10.1
python run_neq_distrib_seed_general.py "/data/chodera/zhangi/perses_benchmark/neq/15/1/" "apo" "$((${LSB_JOBINDEX}-1))" "asn" "lys" 1 # lowercase
