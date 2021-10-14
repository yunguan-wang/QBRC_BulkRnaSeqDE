#python remove_empty_reads_pair.py -RNAfile /project/CRI/Zhu_lab/s171162/Polyploid/seq/rna/YHL_Sample11_R1.fastq.gz /project/CRI/Zhu_lab/s171162/Polyploid/seq/rna/YHL_Sample11_R1.fastq.gz  -out /project/CRI/Zhu_lab/s171162/Polyploid/result/featurecount/YHL_Sample11/
import os
import sys

args=sys.argv
R1seq=args[args.index('-RNAfile')+1]
R2seq=args[args.index('-RNAfile')+2]
output=args[args.index('-out')+1]
def remove_empty_reads(file1,file2,w1,w2):
    with open(file1) as sf1, open(file2) as sf2:
        while True:
            line11=sf1.readline()
            line12=sf1.readline()
            line13=sf1.readline()
            line14=sf1.readline()
            l1=line11+line12+line13+line14
            line21=sf2.readline()
            line22=sf2.readline()
            line23=sf2.readline()
            line24=sf2.readline()
            l2=line21+line22+line23+line24
            if (
                    '' not in l1.split('\n')[0:4] and 
                    '' not in l2.split('\n')[0:4] and 
                    l1.split('\n')[0].replace('/',' ').split(' ')[0] == l2.split('\n')[0].replace('/',' ').split(' ')[0]
            ):
                w1.write(l1)
                w2.write(l2)
            l1=''
            l2=''
            if not line11:
                    break
    sf1.close()
    sf2.close()
                
#convert gzip files to uncompressed files
cmd1='zcat '+R1seq+' > '+output+'/fastq1.org.fq'
cmd2='zcat '+R2seq+' > '+output+'/fastq2.org.fq'
os.system(cmd1)
os.system(cmd2)
#remove empty reads                                                                                                                                          
w1=open(output+'/fastq1.fq','a')
w2=open(output+'/fastq2.fq','a')
remove_empty_reads(output+'/fastq1.org.fq',output+'/fastq2.org.fq',w1,w2)
w1.close()
w2.close()
os.system('rm '+output+'/fastq1.org.fq')
os.system('rm '+output+'/fastq2.org.fq')

