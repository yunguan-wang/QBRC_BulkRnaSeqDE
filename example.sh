#!/bin/bash
#BSUB -W 96:00
#BSUB -n 1
#BSUB -M 128000
#BSUB -e BSUB_output_%j
#BSUB -o BSUB_error_%j

source ~/.bashrc
module load STAR/2.7.4
module load fastqc/0.11.7
module load subread/2.0.3
module load R
###########JOBSTART########################

