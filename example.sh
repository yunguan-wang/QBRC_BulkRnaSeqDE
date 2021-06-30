#!/bin/bash
#SBATCH --job-name=example                                           # job name
#SBATCH --partition=256GB                                       # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --ntasks=32                                                  # number of total tasks
#SBATCH --time=10-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j
source ~/.bash_profile

###########JOBSTART########################

raku /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/BepiTBR_fasta.raku \
--fasta0=/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_data_BepiTBR_fasta/peptide_with_full_length.txt \
--length=15 \
--bepipred2=/project/shared/xiao_wang/projects/Bcell_epitope/code/conda_envs/bp2/bin/activate \
--bepipred1=/project/DPDS/Xiao_lab/shared/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
--LBEEP=/project/shared/xiao_wang/software/LBEEP/ \
--MixMHC2pred=/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=NA \
--dir=/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR_fasta \
--thread=20 \
--keep=false

raku /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/BepiTBR_full.raku \
--full0="mtenstsapaakpkrakaskkstdhpkysdmivaaiqaeknragSSRQSIQKYIKSHYKvgenadsqiklsikrlvttgvlkqtkgvgag\
sfrlaksdepkksvafkktkkeikkvatpkkaskpkkaaskaptkkpkatpvkkakkklaatpkkakkpktvkakpvkaskpkkakpvkpkakssakragkkk" \
--length=15 \
--bepipred2=/project/shared/xiao_wang/projects/Bcell_epitope/code/conda_envs/bp2/bin/activate \
--bepipred1=/project/DPDS/Xiao_lab/shared/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
--LBEEP=/project/shared/xiao_wang/software/LBEEP/ \
--MixMHC2pred=/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=NA \
--dir=/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR_full_new \
--thread=20 \
--keep=false

