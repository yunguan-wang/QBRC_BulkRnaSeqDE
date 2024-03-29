###########  nucleus-specific sbatch wrapper for expression calling  #############
# for other job submission system, you should be able to easily change "sbatch" to the appropriate command
# need 256GB nodes
#
# The instructions (job_expression.pl and expression.pl) must be followed exactly!!! 
#
# input format:
# jobs: the batch job design file, it has 4 columns separated by \t, the first two are fastq files, 
#       the third is the output folder, and the last is "PDX" or "human" or "mouse". 
#       Commented lines ("#" at the front) are skipped. The second fastq file can be left in blank ("")
# example: the demo job submission shell script. The header lines before JOBSTART will be taken for job submission
# index: STAR index
# gtf,thread,disambiguate,count,temp,short: follow those in expression.pl
# n: bundle $n somatic calling jobs into one submission

#perl /project/shared/xiao_wang/software/rnaseqDE/job_expression.pl \
#/project/shared/xiao_wang/software/rnaseqDE/example_data/design4expression.txt \
#/project/shared/xiao_wang/software/rnaseqDE/example.sh \
#/project/shared/xiao_wang/data/hg38/STAR \
#/project/shared/xiao_wang/data/hg38/hg38_genes.gtf \
#32 /project/shared/xiao_wang/software/disambiguate_pipeline \
#count delete N 3

#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($jobs,$example,$index,$gtf,$thread,$disambiguate,$count,$temp,$short,$n)=@ARGV;
my ($line,$line1,@items,$i,$job);

my $path=abs_path($0);
$path=~s/job_expression\.pl//;

open(JOB,$jobs) or die "Cannot find the design file!\n";
$i=0;

while ($line=<JOB>)
{
  $line=~s/(\r|\n)//g;
  if ($line eq "" || $line=~/^#/) {next;}

  if ($i++ % $n==0)
  {
    $job="expression_".$i.".sh";
    open(SCRIPT,">".$job) or die "Cannot write to the shell submission script!\n";

    # write header
    open(HEADER,$example) or die "Cannot find the example shell script!\n";
    while ($line1=<HEADER>)
    {
      if ($line1=~/JOBSTART/) {last;}
      print SCRIPT $line1;
    }
    close(HEADER);
  }

  # write submission job
  @items=split("\t",$line);
  system("rm -f -r ".$items[2]);
  system("mkdir ".$items[2]);
  print SCRIPT "perl ".$path."/expression.pl ".$index.":".$items[0].",".$items[1]." ".$gtf." ".$items[2].
    " ".$thread." ".$items[3]." ".$disambiguate." ".$count." ".$temp." ".$short.
    " 2> ".$items[2]."/job_error.txt 1> ".$items[2]."/job_output.txt\n";

  if ($i % $n==0)
  {
    close(SCRIPT); 
    system("sbatch ".$job);
    unlink($job);
    sleep(1);
  }
}

close(JOB);

if ($i % $n!=0)
{
  close(SCRIPT);
  system("sbatch ".$job);
  unlink($job);
}

