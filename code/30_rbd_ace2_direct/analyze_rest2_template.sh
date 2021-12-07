#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 20:00
#
# Set output file
#BSUB -o 14.%I.analysis.out
#
# Set error file
#BSUB -eo 14.%I.analysis.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu ld-gpu lv-gpu lu-gpu lg-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "14.1[1-4]"

source ~/.bashrc
conda activate perses-rbd-ace2-direct3

old_aa="ASN"
new_aa="LYS"
resid="439"
new=1
rest_dir=/data/chodera/zhangi/perses_benchmark/neq/14/$new
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/15/$new

# Run and analyze rest
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/30_rbd_ace2_direct/
python analyze_cache_general_no_decor.py $rest_dir/ $resid $old_aa $new_aa "${LSB_JOBINDEX}" --length 2

# Make neq directories and copy caches there
mkdir $neq_dir
cp $rest_dir/*.npy $neq_dir
cp $rest_dir/*_apo.pickle $neq_dir
cp $rest_dir/*_complex.pickle $neq_dir

