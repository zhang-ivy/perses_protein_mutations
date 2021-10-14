#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:30 
#
# Set output file
#BSUB -o 14.8.out
#
# Set error file
#BSUB -eo 14.8.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu lx-gpu ly-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=7]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "14.8"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2

old=16
new=9
resid=484
new_aa="GLU"

old_path="/data/chodera/zhangi/perses_benchmark/neq/14/$old/debug/"
#old_path="/data/chodera/zhangi/perses_benchmark/repex/31/$old/0/debug/"
#new_path="/data/chodera/zhangi/perses_benchmark/repex/31/$new/0/"
new_path="/data/chodera/zhangi/perses_benchmark/neq/14/$new/"
rbd_file="4_rbd_ace2_mutant_rbd_tleap_final.pdb"
ace2_file="4_rbd_ace2_mutant_ace2_tleap_final.pdb"

cp "${old_path}${rbd_file}" "${new_path}"0_rbd.pdb
cp "${old_path}${ace2_file}" "${new_path}"0_ace2.pdb

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/21_rbd_ace2_prep/generate_htf.py $new_path $resid $new_aa --rbd_file "${new_path}"0_rbd.pdb --ace2_file "${new_path}"0_ace2.pdb
