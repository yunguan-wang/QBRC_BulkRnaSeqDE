#python remove_empty_reads_single.py -RNAfile /project/CRI/Zhu_lab/s171162/Polyploid/seq/rna/ -out /project/CRI/Zhu_lab/s171162/Polyploid/result/featurecount/YHL_Sample11/ 
import os
import sys
import gzip
args=sys.argv
output=args[args.index('-out')+1]
R1seq=args[args.index('-RNAfile')+1]
def remove_empty_reads(file1,w1):
        with open(file1) as sf1:
                while True:
                        line11=sf1.readline()
                        line12=sf1.readline()
                        line13=sf1.readline()
                        line14=sf1.readline()
                        l1=line11+line12+line13+line14
                        if '' not in l1.split('\n')[0:4]:
                                w1.write(l1)
                        l1=''
                        if not line11:
                                break
        sf1.close()
#convert gzip files to uncompressed files
cmd1='zcat '+R1seq+' > '+output+'/fastq1.org.fq'
os.system(cmd1)
#remove empty reads                                                                                                                                          
w1=open(output+'/fastq1.fq','a')
remove_empty_reads(output+'/fastq1.org.fq',w1)
w1.close()
os.system('rm '+output+'/fastq1.org.fq')
