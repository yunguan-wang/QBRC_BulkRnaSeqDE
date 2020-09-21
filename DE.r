# Differential analysis with DESeq2
# use: R>=3.6

# positional arguments:
#   counts_fn             Input counts table, txt or csv. Genes should be in raw counts.

#   design_fn             Input design table, must be a txt file. First column is condition, second column in contrast. Different contrasts are seperated by ","

# optional arguments:
#   -h, --help            show this help message and exit

#   --output [OUTPUT], -o [OUTPUT]
#     Output path. Default = "./results"

#   --gsea_ref [GSEA_REF], -r [GSEA_REF] 
#     GSEA library path. Must specify full path. Default = "/project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt".

#   --gsea_plot [GSEA_PLOT], -p [GSEA_PLOT] 
#     If GSEA plots will be made. By default no plots are made. Revert with "T"

#   --fccutoff [FCCUTOFF], -f [FCCUTOFF]
#     Log fold change cutoff for volcano plot and the heatmap plot

#   --pcutoff [PCUTOFF], -p [PCUTOFF]
#     Adjusted p-value cutoff for volcano plot.

# Try this:
# Rscript DE.r ./example_data/example_expression.txt ./example_data/example_group.txt -o ./results 
#   -r /project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt -p T

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("EnhancedVolcano", "DESeq2","fgsea")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(packages, rownames(installed.packages())))
}

packages <- c("ggplot2", "dplyr",'tibble','argparse')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

suppressMessages(library('DESeq2'))
suppressMessages(library('ggplot2'))
suppressMessages(library('EnhancedVolcano'))
suppressMessages(library('fgsea'))
suppressMessages(library('dplyr'))
suppressMessages(library('tibble'))
suppressMessages(library('pheatmap'))
suppressMessages(library('argparse'))

# testing
# design = 'Y:/software/rnaseqDE/example_data/example_group.txt'
# cts = 'Y:/software/rnaseqDE/example_data/example_expression.txt'
# output_path = 'Y:/software/rnaseqDE/example_results/'
# gsea_ref = "W:/ref/msigdb/h.all.v7.0.symbols.gmt"

parser <- ArgumentParser(
  description='Differential analysis with DESeq2')

parser$add_argument(
  'counts_fn', type="character",
  help='Input counts table, txt or csv. Genes should be in raw counts.')

parser$add_argument(
  'design_fn', type="character",
  help='Input design table, must be a txt file. First column is condition, second column in contrast. Different contrasts are seperated by ","')

parser$add_argument(
  '--output','-o',nargs='?', default='./results',
  help='Output path. Default = "./results"')

parser$add_argument(
  '--gsea_ref','-r',nargs='?', 
  default='/project/shared/brugarolas_wang_xiao/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt',
  help='GSEA library path. Default = "/project/shared/brugarolas_wang_xiao/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt".')

parser$add_argument(
  '--gsea_plot','-g',nargs='?', default='N',
  help='If GSEA plots will be made. By default no plots are made. Revert with "T"')

parser$add_argument(
  '--fccutoff','-f',nargs='?', default=2, type='double',
  help='Log fold change cutoff for volcano plot.')

parser$add_argument(
  '--pcutoff','-p',nargs='?', default=1e-3, type='double',
  help='Adjusted p-value cutoff for volcano plot.')

args <- parser$parse_args()
cts <- args$counts_fn
design <- args$design_fn
output_path <- args$output
gsea_ref <- args$gsea_ref
gsea_plot <- args$gsea_plot
fc_cutoff <- args$fccutoff
pval_cutoff <- args$pcutoff

# --------
# Preprocessing input data

# gene names must not have repeats
design <- read.table(design,stringsAsFactors = F,header=T, sep='\t',row.names = 1)
cts <- read.table(cts, stringsAsFactors = F,header=T, sep='\t', row.names = 1)

# Make output file and set path to it
dir.create(file.path(output_path), showWarnings = FALSE)
setwd(file.path(output_path))

# Setting up contrasts
contrast_groups <- unique(design)
analysis = list()
j=1
for (i in 1:dim(contrast_groups)[1]) {
  t <- contrast_groups[i,1]
  refs <- strsplit(contrast_groups[i,2],',')[[1]]
  for (ref in refs) {
    analysis[[j]] <- list(t,ref)
    j <- j + 1
  }
} 

# Align dataset with design
design <- subset(design,row.names(design) %in% colnames(cts))
cts <- cts[,row.names(design)]
colnames(design) <- c('condition','Contrasts')

# --------
# DEseq2 analysis
dds <- DESeqDataSetFromMatrix(
  countData = cts, colData = design, design = ~ condition)
# Simply checking for not expressed genes, filter at sum of counts of all
# samples >= 10
keep <- rowSums(counts(dds)) >= 10
dds = dds[keep,]

# PCA based on DESeq2 normalized counts
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(
  pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave('PCA_plot.pdf',height = 5, width = 5)

# Distance plot
corr <- cor(assay(vsd), method='pearson')
pheatmap(
  corr,annotation_col = design['condition'],show_rownames=F,cluster_cols = T,
  filename = "Sample distance heatmap.pdf",sep = '')

# Write normalized count to results.
export_counts <- assay(vsd)
write.table(export_counts,'deseq2 vst counts.csv', sep=',')

# Looping through all contrasts
for (c in analysis){
  target <- c[[1]]
  ref <- c[[2]]
  # Todo: Implement grouping joining feature.
  # This should allow comparing group A with B+C combined.
  if (grepl('|', ref, fixed = TRUE)) {
    print("Pooling conditions not supported!")
    next
  }
  
  # Comparison is defined as target vs ref
  res = results(dds,contrast = c("condition",target,ref))
  output_prefix = paste(target,"_vs_",ref, '/', sep = '')
  dir.create(file.path(output_prefix), showWarnings = FALSE)
  # writing results
  write.table(as.data.frame(res),
            file=paste(output_prefix,'/DEG.txt',sep=""),sep='\t') # by wtwt5237
  
  # Volcano plot
  EnhancedVolcano(
    res, lab = rownames(res),x = 'log2FoldChange',y = 'pvalue',xlim = c(-5, 5),
    ylim = c(0,10), pCutoff = pval_cutoff,title = '', FCcutoff = fc_cutoff,
    subtitle = output_prefix
    )
  ggsave(paste(output_prefix,"/DEG_Volcano_plot.pdf",sep = '')) # by wtwt5237
  
  # Heatmap
  res_heatmap = res[res$padj<=pval_cutoff & !is.na(res$padj),] # by wtwt5237
  # by wtwt5237 -start
  tops_abs = rownames(res_heatmap[order(abs(res_heatmap$stat), decreasing = T),])[1:min(200,dim(res_heatmap)[1])]
  if (length(tops_abs) == 0) {
    next
  }
  heatmap_mats = export_counts[tops_abs,]

  # by wtwt5237 - end
  pheatmap(
    heatmap_mats,annotation_col = design['condition'],show_rownames=T, # by wtwt5237
    scale = 'row',cluster_cols = F, # by wtwt5237 
    filename = paste(output_prefix,"/DEG_heatmap.pdf",sep = '')) # by wtwt5237
  
  # ========
  # GSEA analysis
  # pre-processing Deseq results
  res$gene = row.names(res)
  res2 = res[,c("gene","stat")]
  ranks <- deframe(res2)
  pathway_kegg = gmtPathways(gsea_ref)
  
  # GSEA PreRank using LFC stat in DESeq results, which is basically logFC.
  fgseaRes_kegg <- fgseaMultilevel(pathways=pathway_kegg, stats=ranks)
  fgseaResTidy_c <- fgseaRes_kegg %>%
    as_tibble() %>%
    arrange(desc(NES))

  # save pathway data
  df_pathway <- as.data.frame(fgseaResTidy_c)
  df_pathway$leadingEdge <- unlist(lapply(df_pathway$leadingEdge, function (x){
    paste(x,collapse=', ')}))
  write.table(df_pathway, file=paste(output_prefix,'/pathways.txt',sep=""), sep='\t') # by wtwt5237
  
  # Todo: Added a handle to turn this off in argument.
  # Bar plot of all pathways passing padj <=0.25 mark. 
  fgseaResTidy_c = fgseaResTidy_c[fgseaResTidy_c$padj<=0.25,]
  ggplot(
    fgseaResTidy_c, aes(reorder(pathway, NES), NES)
  ) + geom_col(aes(fill=padj<0.25)) + coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
          title="KEGG from GSEA") +
    theme_minimal()
  ggsave(paste(output_prefix, "/GSEA.pdf",sep=""), width=9,height=16) # by wtwt5237
  
  # Make all GSEA plots
  if (gsea_plot == 'T'){
    if(! dir.exists('GSEA_plots')){
      dir.create('GSEA_plots')
    }
    
    if (!dir.exists(paste(output_prefix,"/GSEA_plots",sep=""))) # by wtwt5237
    dir.create(paste(output_prefix,"/GSEA_plots",sep="")) # by wtwst5237

    for (pathway_name in fgseaResTidy_c$pathway){
      plotEnrichment(pathway_kegg[[pathway_name]], ranks)
      ggsave(paste(output_prefix,'/GSEA_plots/', pathway_name, '.pdf',sep='')) # by wtwt5237
    }
  }
}

warnings() # by wtwt5237
