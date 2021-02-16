#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o  13.63.htf.%I.out
#
# Set error file
#BSUB -eo 13.63.htf.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu lt-gpu lu-gpu lv-gpu lw-gpu lx-gpu ly-gpu" 
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=2]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "13.63[1]"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
python generate_htfs.py "/data/chodera/zhangi/perses_benchmark/neq/13/63/" "87" "ALA" "ARG" "${LSB_JOBINDEX}" "../../input/mmc2_barnase_R87A.pdb" "../../input/mmc2_barstar.pdb"

