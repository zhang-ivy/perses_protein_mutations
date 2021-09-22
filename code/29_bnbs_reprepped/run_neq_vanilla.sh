#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 6:00 
#
# Set output file
#BSUB -o  23.1.complex.%I.out
#
# Set error file
#BSUB -eo 23.1.complex.%I.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu ld-gpu lv-gpu lu-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "23.1.com[1-100]"
# wait for another job to be completed

old_aa="ASN"
new_aa="LYS"
resid="439"
new=1
htf_dir=/data/chodera/zhangi/perses_benchmark/neq/22/$new
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/23/$new/
phase='complex'
eq_length=1
neq_length=1

source ~/.bashrc
module load cuda/10.2
conda activate perses-nonstandard-aa

htf_file="$new"_"$phase".pickle
if [ ! -f $neq_dir"$htf_file" ]; then
    cp $htf_dir/$htf_file $neq_dir
fi

cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/29_bnbs_reprepped/
python run_neq_distrib_seed_general_vanilla.py $neq_dir $phase "$((${LSB_JOBINDEX}-1))" $old_aa $new_aa $eq_length $neq_length # lowercase
