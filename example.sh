#!/bin/bash
#BSUB -W 24:00
#BSUB -n 4
#BSUB -M 128000
#BSUB -e BSUB_output_%j
#BSUB -o BSUB_error_%j

source ~/.bashrc
conda activate p3s
module load STAR/2.7.4
module load fastqc/0.11.7
module load subread/2.0.3

###########JOBSTART########################

