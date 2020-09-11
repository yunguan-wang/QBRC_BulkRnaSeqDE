# prerequisite in path: python, featureCounts, STAR (>=2.7.2b), fastqc, disambiguate (use conda env by Yunguan)
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

# (1) STAR index:fastq1,fastq2. fastq2 is optional. Fastq files must be gzipped
# (2) gtf file
# (3) output folder
# (4) number of threads to use
# (5) pdx or not. "PDX" or "human" or "mouse". if "PDX", can only handle paired-end sequencing reads
# (6) disambiguate path: prepared by Yunguan (yunguan.wang@utsouthwestern.edu), 
#     default: /project/shared/xiao_wang/software/disambiguate_pipeline
# (7) count: "rpkm" or "count" (integer)? This must match with the input for summarize4expression.py
#
#perl /project/bioinformatics/Xiao_lab/shared/neoantigen/code/expression/expression.pl \
#/project/shared/xiao_wang/data/hg38/STAR:\
#/archive/BICF/shared/Kidney/rna/RAW/SAM24297102_1_R1.fastq.gz,\
#/archive/BICF/shared/Kidney/rna/RAW/SAM24297102_1_R2.fastq.gz \
#/project/shared/xiao_wang/data/hg38/hg38_genes.gtf \
#/project/SCCC/Wang_lab/shared/tmp \
#32 human \
#/project/shared/xiao_wang/software/disambiguate_pipeline
#count

my ($bam_file,$gtf,$output_folder,$thread,$pdx,$disambiguate,$count)=@ARGV;
my ($i,$dir,$path,$index,$fastq,$fastq_new);
my $parameters=" --primary -O -t exon -g transcript_id -T ".$thread." --largestOverlap --minOverlap 3 --ignoreDup -p -P -B -C";

#################  clean data  ###############################

system("rm -f -r ".$output_folder);
system("mkdir ".$output_folder);

$path=abs_path($0);
$path=~s/expression\.pl//;

$bam_file=~/^(.*)\:(.*)$/;
$index=$1;
$fastq=$2;
$fastq=~s/,/ /;

#############33  preprocess fastq files  ####################

$i=1;
$fastq_new="";

foreach (split(" ",$fastq))
{
  system("zcat ".$_." > ".$output_folder."/zcat_".$i.".fastq");
  $fastq_new.=$output_folder."/zcat_".$i++.".fastq ";
}
$fastq=$fastq_new;

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
  unlink($output_folder."/zcat_1.fastq");
  unlink($output_folder."/zcat_2.fastq");
}

# STAR alignment
system("STAR --runThreadN ".$thread." --genomeDir ".$index." --readFilesIn ".$fastq." --outFileNamePrefix ".$output_folder.
  "/ --outSAMtype BAM Unsorted");
unlink($output_folder."/Unmapped.out.mate1");
unlink($output_folder."/Unmapped.out.mate2");

# featureCounts
$bam_file=$output_folder."/Aligned.out.bam";
unless (-f $gtf) {die "Error: Gtf annotation file doesn't exists!\n";}
unless (-f $bam_file) {die "Error: RNA-Seq bam file doesn't exists!\n";}

system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts ".$bam_file." ".$parameters." -s 0");
system("rm -f ".$output_folder."/aligned.sam");
unless (-f $output_folder."/transcript.featureCounts") {die "Error: featureCounts failed!\n";}
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts_stranded ".$bam_file." ".$parameters." -s 1");
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts_rev_stranded ".$bam_file." ".$parameters." -s 2");

# rpkm
system("Rscript ".$path."/script/rpkm.R ".$output_folder."/transcript.featureCounts ".$count);

##############  cleanup  #######################

unlink($output_folder."/Log.out");
unlink($output_folder."/Log.progress.out");
unlink($output_folder."/SJ.out.tab");
unlink($output_folder."/*summary");
unlink($bam_file);
unlink($output_folder."/zcat_1.fastq");
unlink($output_folder."/zcat_2.fastq");
unlink($output_folder."/fastq1.fastq.human");
unlink($output_folder."/fastq2.fastq.human");

