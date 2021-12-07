#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 25:00 
#
# Set output file
#BSUB -o  complex.out
#
# Set error file
#BSUB -eo complex.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lw-gpu lu-gpu lg-gpu lt-gpu"
#BSUB -q gpuscatter
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=7]"
#BSUB -n 2 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "32.0.0"

source ~/.bashrc
module load cuda/10.2
OMPI_MCA_opal_cuda_support=true
conda activate perses-rbd-ace2-direct2

new=0
outdir=/data/chodera/zhangi/perses_benchmark/repex/32/0/$new/
phase='complex'
lambda_start=0.0
lambda_end=0.1
n_states=12
n_cycles=5000

build_mpirun_configfile --configfilepath configfile_complex --hostfilepath hostfile_complex "python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/32_rbd_ace2_direct_repex/run_h_repex_chunk.py $outdir $phase $lambda_start $lambda_end $n_states $n_cycles"
mpiexec.hydra -f hostfile_complex -configfile configfile_complex
