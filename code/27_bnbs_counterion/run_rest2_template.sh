#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 12:00
#
# Set output file
#BSUB -o 17.%I.rest.out
#
# Set error file
#BSUB -eo 17.%I.rest.stderr
#
# Specify node group
#BSUB -m "lx-gpu ly-gpu lw-gpu lu-gpu lv-gpu ld-gpu lp-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "17.1r[1-4]"

source ~/.bashrc
module load cuda/10.2
conda activate perses-counterion-2

old_aa="ASN"
new_aa="LYS"
resid="439"
new=1
rest_dir=/data/chodera/zhangi/perses_benchmark/neq/17/$new
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/18/$new

# Run and analyze rest
python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/27_bnbs_counterion/generate_rest2_cache_auto_lang_general_1ps_custom.py $rest_dir/ $resid $old_aa $new_aa "${LSB_JOBINDEX}" 2000 0.4 24


