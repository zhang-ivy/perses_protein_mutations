#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:30 
#
# Set output file
#BSUB -o  17.43.%I.out
#
# Set error file
#BSUB -eo 17.43.%I.stderr
#
# Specify node group
#BSUB -m "lx-gpu ly-gpu lw-gpu lu-gpu lv-gpu ld-gpu lp-gpu lg-gpu lt-gpu" 
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=4]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "17.43"

source ~/.bashrc
module load cuda/10.2
conda activate perses-counterion-2

project=17
dir=43
resid=38
new_aa=TRP
protein_path="/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/mmc2_barstar.pdb"
ligand_path="/home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/input/mmc2_barnase.pdb"

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/27_bnbs_counterion/generate_htfs.py /data/chodera/zhangi/perses_benchmark/neq/$project/$dir/ $resid $new_aa --input_file $protein_path --ligand_file $ligand_path

