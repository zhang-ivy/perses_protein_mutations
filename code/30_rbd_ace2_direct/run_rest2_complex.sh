#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 20:00
#
# Set output file
#BSUB -o 14.113.%I.rest.out
#
# Set error file
#BSUB -eo 14.113.%I.rest.stderr
#
# Specify node group
#BSUB -q gpuscatter
#BSUB -m "ly-gpu lx-gpu lw-gpu lu-gpu lg-gpu"
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=10]"
#BSUB -n 4 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "14.113r[3-4]"

source ~/.bashrc
module load cuda/10.2
conda activate perses-rbd-ace2-direct

old_aa=ASN
new_aa=TYR
resid=501
new=113
rest_dir=/data/chodera/zhangi/perses_benchmark/neq/14/$new
neq_dir=/data/chodera/zhangi/perses_benchmark/neq/15/$new

# Run and analyze rest
build_mpirun_configfile --configfilepath configfile_complex_"${LSB_JOBINDEX}" --hostfilepath hostfile_complex_"${LSB_JOBINDEX}" "python $rest_dir/generate_rest2_cache_auto_lang_general_1ps_custom.py $rest_dir/ $resid $old_aa $new_aa ${LSB_JOBINDEX} 1200 0.5"
mpiexec.hydra -f hostfile_complex_"${LSB_JOBINDEX}" -configfile configfile_complex_"${LSB_JOBINDEX}"
