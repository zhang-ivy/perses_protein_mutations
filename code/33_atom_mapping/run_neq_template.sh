#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  15.303.complex.%I.out
#
# Set error file
#BSUB -eo 15.303.complex.%I.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu ld-gpu lv-gpu lu-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "15.303.c[1-100]"

new=303
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/15/$new/
phase=complex
eq_length=1
neq_length=1
chain_A=0
chain_B=2

source ~/.bashrc
conda activate openmm-dev

cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/31_rest_over_protocol/
python run_neq_vanilla.py $neq_dir $phase "$((${LSB_JOBINDEX}-1))" $eq_length $neq_length $chain_A $chain_B
