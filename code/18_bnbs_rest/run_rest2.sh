#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 8:00 
#
# Set output file
#BSUB -o  12.52.4.out
#
# Set error file
#BSUB -eo 12.52.4.stderr
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
#BSUB -J "12.52.4"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
#python generate_rest2_cache_bnbs_auto_lang.py "/data/chodera/zhangi/perses_benchmark/neq/12/39/" "apo" "ALA" 0 5 1200 "forward" 3 4.0 0.2 12
python generate_rest2_cache_bnbs_auto_lang_general.py "/data/chodera/zhangi/perses_benchmark/neq/12/52/" "76" "GLU" "ALA" 4 "../../input/mmc2_barstar.pdb"
