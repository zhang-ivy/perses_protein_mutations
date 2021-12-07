#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 5:59
#
# Set output file
#BSUB -o 14.113.%I.rest.out
#
# Set error file
#BSUB -eo 14.113.%I.rest.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu ld-gpu lv-gpu lu-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=7]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "14.113r[1-2]"

source ~/.bashrc
conda activate perses-rbd-ace2-direct3

old_aa=ASN
new_aa=TYR
resid=501
new=113
rest_dir=/data/chodera/zhangi/perses_benchmark/neq/14/$new
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/15/$new

# Run and analyze rest
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/30_rbd_ac32_direct/
python generate_rest2_cache_auto_lang_general_1ps_custom.py $rest_dir/ $resid $old_aa $new_aa "${LSB_JOBINDEX}" 1200 0.3 12


