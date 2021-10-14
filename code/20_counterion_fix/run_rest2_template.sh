#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 20:00
#
# Set output file
#BSUB -o 16.4.%I.out
#
# Set error file
#BSUB -eo 16.4.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lu-gpu lp-gpu lg-gpu lt-gpu lv-gpu lw-gpu lx-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "16.4[1-2]"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2

old_aa="ASN"
new_aa="LYS"
resid="439"
new=4
rest_dir=/data/chodera/zhangi/perses_benchmark/neq/16/$new

# Run and analyze rest
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/20_counterion_fix/
python generate_rest2_cache_auto_lang_general.py $rest_dir/ $resid $old_aa $new_aa "${LSB_JOBINDEX}"
python analyze_cache_general.py $rest_dir/ $resid $old_aa $new_aa "${LSB_JOBINDEX}"


