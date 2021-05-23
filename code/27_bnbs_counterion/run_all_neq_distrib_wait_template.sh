#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  18.1.complex.%I.out
#
# Set error file
#BSUB -eo 18.1.complex.%I.stderr
#
# Specify node group
#BSUB -m "lx-gpu ly-gpu lw-gpu lu-gpu lv-gpu ld-gpu lp-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "18.1.com[1-100]"
# wait for another job to be completed
#BSUB -w "done(79112)"

old_aa="ASN"
new_aa="LYS"
resid="439"
new=1
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/18/$new/
phase='complex'
length=1

source ~/.bashrc
module load cuda/10.2
conda activate perses-counterion-2

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/27_bnbs_counterion/run_neq_distrib_seed_general.py $neq_dir $phase "$((${LSB_JOBINDEX}-1))" $old_aa $new_aa $length # lowercase
