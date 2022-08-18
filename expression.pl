# prerequisite in path: python (v2 or v3), featureCounts, R, STAR (>=2.7.2b), fastqc, disambiguate (use conda env by Yunguan)
# can now handle hg38, mm10, sacCer3
# need 256 GB nodes
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

# (1) STAR index:fastq1,fastq2. fastq2 is optional. Fastq files must be gzipped
# (2) gtf file
# (3) output folder
# (4) number of threads to use
# (5) pdx or not. "PDX" or "human" or "mouse" or "yeast". if "PDX", can only handle paired-end sequencing reads
# (6) disambiguate path: prepared by Yunguan (yunguan.wang@utsouthwestern.edu), 
#     default: /project/shared/xiao_wang/software/disambiguate_pipeline
# (7) count: "rpkm" or "count" (unnormalized and integer valued)? For running DE.r, use "count"
# (8) temp: "keep" or "delete" temporary files (default is "delete")
# (9) short: "Y" or "N". Default is "N"
#
#perl /project/shared/xiao_wang/software/rnaseqDE/expression.pl \
#/project/shared/xiao_wang/data/hg38/STAR:\
#/archive/BICF/shared/Kidney/rna/RAW/SAM24297102_1_R1.fastq.gz,\
#/archive/BICF/shared/Kidney/rna/RAW/SAM24297102_1_R2.fastq.gz \
#/project/shared/xiao_wang/data/hg38/hg38_genes.gtf \
#/project/DPDS/Wang_lab/shared/tmp \
#32 human \
#/project/shared/xiao_wang/software/disambiguate_pipeline \
#count delete N

my ($bam_file,$gtf,$output_folder,$thread,$pdx,$disambiguate,$count,$temp,$short)=@ARGV;
my ($short_parameters,$fc_parameters,$i,$dir,$path,$index,$fastq,$fastq_new);

#################  set path, clean data, and set parameters  ###############################

$path=abs_path($0);
$path=~s/expression\.pl//;

$bam_file=~/^(.*)\:(.*)$/;
$index=$1;
$fastq=$2;
$fastq=~s/,/ /;

# delete emptoy reads
if (! -d $output_folder) {system("mkdir ".$output_folder);}
system("mkdir ".$output_folder."/clean_fastq");

if ($fastq!~/\s$/)
{
  system("python ".$path."/script/remove_empty_reads_pair.py ".
    "-out ".$output_folder."/clean_fastq -RNAfile ".$fastq);
  $fastq=$output_folder."/clean_fastq/fastq1.fq ".$output_folder."/clean_fastq/fastq2.fq";
  $fc_parameters=" --primary -O -t exon -g transcript_id -T ".$thread." --largestOverlap --minOverlap 3 --ignoreDup -p -P -B -C";
}else
{
  system("python ".$path."/script/remove_empty_reads_single.py ".
    "-out ".$output_folder."/clean_fastq -RNAfile ".$fastq);
  $fastq=$output_folder."/clean_fastq/fastq1.fq";  
  $fc_parameters=" --primary -O -t exon -g transcript_id -T ".$thread." --largestOverlap --minOverlap 3 --ignoreDup";
}

# set star and fc parameters
if ($short eq "N") 
{
  $short_parameters="";
}else
{
  $short_parameters="--seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 ".
    "--outFilterMatchNmin 10 --outFilterMultimapNmax 10 ";
  $fc_parameters.=" -M "; 
}

#################  quality check  ############################

system("mkdir ".$output_folder."/fastqc");
system("fastqc -o ".$output_folder."/fastqc --extract -t ".$thread." -q -d ".$output_folder."/fastqc ".$fastq);

$i=1;
opendir(DIR,$output_folder."/fastqc");
foreach $dir (readdir(DIR))
{
  if ($dir!~/fastqc$/) {next;}
  system("mv ".$output_folder."/fastqc/".$dir."/summary.txt ".$output_folder."/qc_summary_".$i++.".txt");
}
close(DIR);
system("rm -f -r ".$output_folder."/fastqc");

#################  expression analyses  ######################

# PDX model
if ($pdx eq "PDX")
{
  system("source activate ".$disambiguate."/conda_env;".
    "python ".$disambiguate."/ngs_disambiguate.py -o ".$output_folder.
    " -i ".$output_folder."/disambiguate -a star -r \"".$index."|".$index."/mouse\"".
    " -n ".$thread." -b ".$path."/script/bam2fastq.pl ".$fastq.";source deactivate");
  system("rm -f -r ".$output_folder."/alignment_human");
  system("rm -f -r ".$output_folder."/alignment_mouse");
  system("rm -f -r ".$output_folder."/fastq1.fastq.mouse");
  system("rm -f -r ".$output_folder."/fastq2.fastq.mouse");
  $fastq=$output_folder."/fastq1.fastq.human ".$output_folder."/fastq2.fastq.human";
  system("rm -f -r ".$output_folder."/clean_fastq");
}

# STAR alignment
system("STAR --runThreadN ".$thread." --genomeDir ".$index." --readFilesIn ".$fastq." --quantMode TranscriptomeSAM GeneCounts ".
  "--outFileNamePrefix ".$output_folder."/ --outSAMtype BAM Unsorted ".$short_parameters);
system("mv ".$output_folder."/ReadsPerGene.out.tab ".$output_folder."/STAR_gene_counts.txt");
unlink($output_folder."/Unmapped.out.mate1");
unlink($output_folder."/Unmapped.out.mate2");
system("rm -f -r ".$output_folder."/clean_fastq");
unlink($output_folder."/fastq1.fastq.human");
unlink($output_folder."/fastq2.fastq.human");

# featureCounts
$bam_file=$output_folder."/Aligned.out.bam";
unless (-f $gtf) {die "Error: Gtf annotation file doesn't exists!\n";}
unless (-f $bam_file) {die "Error: RNA-Seq bam file doesn't exists!\n";}

system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts ".$bam_file." ".$fc_parameters." -s 0");
system("rm -f ".$output_folder."/aligned.sam");
unless (-f $output_folder."/transcript.featureCounts") {die "Error: featureCounts failed!\n";}
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts_stranded ".$bam_file." ".$fc_parameters." -s 1");
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts_rev_stranded ".$bam_file." ".$fc_parameters." -s 2");

# rpkm
system("Rscript ".$path."/script/rpkm.R ".$output_folder."/transcript.featureCounts ".$count);

##############  cleanup  #######################

if ($temp eq "delete")
{
  unlink($output_folder."/Log.out");
  unlink($output_folder."/Log.progress.out");
  unlink($output_folder."/SJ.out.tab");
  unlink($output_folder."/*summary");
  unlink($bam_file);
}

