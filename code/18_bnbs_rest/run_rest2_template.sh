#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 12:00
#
# Set output file
#BSUB -o 12.%I.out
#
# Set error file
#BSUB -eo 12.%I.stderr
#
# Specify node group
#BSUB -m "ld-gpu lg-gpu lp-gpu lt-gpu lu-gpu lv-gpu lw-gpu lx-gpu ly-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "12.99[1-4]"

source ~/.bashrc
module load cuda/10.1
conda activate perses-sims-oct2020

old=58
new=99
old_aa="LYS"
new_aa="ALA"
resid="27"

# Make dirs and copy r-htfs for rest
cd /data/chodera/zhangi/perses_benchmark/neq/12/
cp $old/$old\_apo_0.pickle $new/$new\_apo_0.pickle
cp $old/$old\_apo_1.pickle $new/$new\_apo_1.pickle
cp $old/$old\_complex_0.pickle $new/$new\_complex_0.pickle
cp $old/$old\_complex_1.pickle $new/$new\_complex_1.pickle

# Run and analyze rest
cd /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/18_bnbs_rest
python generate_rest2_cache_bnbs_auto_lang_extreme_3.py "/data/chodera/zhangi/perses_benchmark/neq/12/$new/" $resid $old_aa $new_aa "${LSB_JOBINDEX}" "../../input/mmc2_barnase.pdb"
python analyze_cache_general.py "/data/chodera/zhangi/perses_benchmark/neq/12/$new/" $resid $old_aa $new_aa "${LSB_JOBINDEX}"

# Make neq directories and copy caches there
cd /data/chodera/zhangi/perses_benchmark/neq/13/
mkdir $new
cp /data/chodera/zhangi/perses_benchmark/neq/12/$new/*.npy /data/chodera/zhangi/perses_benchmark/neq/13/$new/
cp $old/$old\_apo.pickle $new/$new\_apo.pickle
cp $old/$old\_complex.pickle $new/$new\_complex.pickle

