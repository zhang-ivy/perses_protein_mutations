#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00 
#
# Set output file
#BSUB -o  complex.out
#
# Set error file
#BSUB -eo complex.stderr
#
# Specify node group
#BSUB -m "lt-gpu lp-gpu lg-gpu lv-gpu ld-gpu lu-gpu lx-gpu ly-gpu lw-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=7]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "3.0.com"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/24_rbd_ace2_repex/
python run_h_repex.py "/data/chodera/zhangi/perses_benchmark/repex/31/3/0/" "complex"

