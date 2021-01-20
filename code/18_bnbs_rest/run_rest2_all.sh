#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 5:00 
#
# Set output file
#BSUB -o  12.43.%I.out
#
# Set error file
#BSUB -eo 12.43.%I.stderr
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
#BSUB -J "12.43[1-4]"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020
python generate_rest2_cache_bnbs_auto_lang_general.py "/data/chodera/zhangi/perses_benchmark/neq/12/43/" "38" "PHE" "TRP" "${LSB_JOBINDEX}" "../../input/mmc2_barstar_W38F.pdb"

