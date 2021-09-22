#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:30 
#
# Set output file
#BSUB -o 22.8.out
#
# Set error file
#BSUB -eo 22.8.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu ld-gpu lv-gpu lu-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=4]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "22.8"

source ~/.bashrc
module load cuda/10.2
conda activate perses-nonstandard-aa

new=9
resid=484
new_aa="GLU"
new_path="/data/chodera/zhangi/perses_benchmark/neq/22/$new/"

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/29_bnbs_reprepped/generate_htfs.py $new_path $resid $new_aa
