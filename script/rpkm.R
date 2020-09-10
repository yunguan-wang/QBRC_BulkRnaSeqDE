file0=commandArgs(trailingOnly=TRUE)[1]
count=commandArgs(trailingOnly=TRUE)[2]

assigned_all=list() # for determining strandedness

rpkm<-function(x,count)
{
  if (sum(x$expression)==0) {stop("Error: gene expression count is 0!")}
  if (count=="count") 
  {
    x$expression=round(x$expression)
    return(x)
  }
  x$expression=x$expression/sum(x$expression)*10e6
  x$expression=x$expression/x$Length*1000
  x
}

for (suffix in c("","_stranded","_rev_stranded"))
{
  # count file
  file=paste(file0,suffix,sep="")
  
  transcript_exp=read.table(file,stringsAsFactors = F,header=T)
  colnames(transcript_exp)[7]="expression"
  transcript_exp=rpkm(transcript_exp,count)
  
  write.table(transcript_exp,file=file,sep="\t",quote=F,row.names = F)
  
  # strandedness
  assigned=read.table(paste(file,".summary",sep=""),header=T,sep="\t")[1,2]
  assigned_all[[paste("dummy_",suffix,sep="")]]=assigned
}

# check strandedness
strand="unknown"
if ((assigned_all$dummy__stranded+assigned_all$dummy__rev_stranded)/assigned_all$dummy_>0.9) {strand="unstranded"}
if (assigned_all$dummy__stranded/assigned_all$dummy__rev_stranded>10) {strand="stranded"}
if (assigned_all$dummy__rev_stranded/assigned_all$dummy__stranded>10) {strand="rev_stranded"}

write.table(paste(strand,count),file=paste(file0,"_parameters",sep=""),
            sep="\t",quote=F,row.names = F,col.names = F)
