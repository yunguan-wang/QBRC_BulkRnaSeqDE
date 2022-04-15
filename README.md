### Dependencies
Reads were aligned to reference (GRCh38) with ‘STAR’ (v2.7.2b). 

Gene counts were quantified with ‘FeatureCounts’ (v1.6.4).

Differential gene expression analysis was performed using the R package ‘DEseq2’(v1.26).

GSEA statistical analysis was carried out with the R package ‘fgsea’ (v1.14.0).


### Updates:

12/08/2020: Added handling of short read alignment
03/02/2021: Added STAR gene counts
06/29/2021: Clarified documentation in a number of places
01/25/2022: Writing stderr and stdout to job folders

In this folder

design4expression.py -> job_expression.pl (expression.pl) -> summarize4expression.py -> DE.r

---------------------------------------

File: example.sh

Description: This shell script is an example job submission script to be used by job_expression.pl

---------------------------------------

File: design4expression.py

Description: This python3 script take samples from manifest file, and create two design files, one for job_expression.pl, another for summarize4expression.py

---------------------------------------

File: job_expression.pl

Description: This perl script is a wrapper around expression.pl for easy submission of the expression analysis jobs

---------------------------------------

File: expression.pl

Description: This perl script executes QC, alignment, and counting for raw RNA-seq data

----------------------------------------

File: summarize4expression.py

Description: This python3 script examines the output from expression.pl for all analyzed samples, and summarizes the expression data into a matrix table 

-----------------------------------------

File: DE.r

Description: This R script takes the output of summarize4expression.py and performs regular RNA differential analyses

-----------------------------------------

Folder: script

Description: This folder stores dependency scripts for the main scripts

-----------------------------------------

Folder: example_data

Description: This folder stores example input data for the pipeline

-----------------------------------------

Folder: results

Description: This folder stores the example output of the pipeline
