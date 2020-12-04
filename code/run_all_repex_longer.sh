#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 50:00 
#
# Set output file
#BSUB -o repex-T42A.complex.out
#
# Set error file
#BSUB -eo repex-T42A.complex.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "repex-T42A_complex[13-15]"

module load cuda/10.1
python run_repex_longer.py "/data/chodera/zhangi/perses_benchmark/repex/"${LSB_JOBINDEX}"/0/" "complex"

