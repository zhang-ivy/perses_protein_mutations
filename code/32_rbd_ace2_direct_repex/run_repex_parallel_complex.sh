#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 50:00 
#
# Set output file
#BSUB -o  complex.out
#
# Set error file
#BSUB -eo complex.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=20]"
#BSUB -n 1 -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "32.0.0"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2-direct2

new=0
outdir=/data/chodera/zhangi/perses_benchmark/repex/32/1/$new/
phase=complex
lambda_start=0.0
lambda_end=0.1
n_states=12
n_cycles=5000

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/32_rbd_ace2_direct_repex/run_h_repex_chunk.py $outdir $phase $lambda_start $lambda_end $n_states $n_cycles
