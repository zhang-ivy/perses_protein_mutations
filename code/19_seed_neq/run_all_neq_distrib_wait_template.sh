#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o  13.101.complex.%I.out
#
# Set error file
#BSUB -eo 13.101.complex.%I.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "13.101.com[1-100]"
# wait for another job to be completed

new=107
old_aa=ASP
new_aa=ALA
phase='complex'
length=1

source ~/.bashrc
conda activate perses-sims-oct2020
module load cuda/10.1
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/19_seed_neq/
python run_neq_distrib_seed_general.py "/data/chodera/zhangi/perses_benchmark/neq/13/$new/" $phase "$((${LSB_JOBINDEX}-1))" $old_aa $new_aa $length # lowercase
