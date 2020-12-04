#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 10:00 
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=10]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "ALA_SER"

module load cuda/10.1
export WT="ALA"
export PROPOSED="SER"
python run_aa.py "$WT" "$PROPOSED" &> "${WT}_${PROPOSED}_solvent.log"
