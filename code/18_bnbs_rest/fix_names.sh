#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00 
#
# Set output file
#BSUB -o  12.1.fix.names.out
#
# Set error file
#BSUB -eo 12.1.fix.names.stderr
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "fix names"
#
# wait for another job to be completed
#BSUB -w "done(18988446)"

source ~/.bashrc
conda activate perses-sims
cd /data/chodera/zhangi/perses_benchmark/neq/12/1/
cp *nc old/
mv 1_apo_ala_5ns_checkpoint.nc 1_apo_thr_5ns_checkpoint.nc
mv 1_apo_ala_5ns.nc 1_apo_thr_5ns.nc
mv 1_complex_ala_5ns_checkpoint.nc 1_complex_thr_5ns_checkpoint.nc
mv 1_complex_ala_5ns.nc 1_complex_thr_5ns.nc

cp old/1_apo_thr_5ns_checkpoint.nc 1_apo_ala_5ns_checkpoint.nc
cp old/1_apo_thr_5ns.nc 1_apo_ala_5ns.nc
cp old/1_complex_thr_5ns_checkpoint.nc 1_complex_ala_5ns_checkpoint.nc
cp old/1_complex_thr_5ns.nc 1_complex_ala_5ns.nc
