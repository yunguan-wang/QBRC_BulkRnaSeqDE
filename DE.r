# Differential analysis with DESeq2 >= 1.26
# use: R>=3.6

# positional arguments:
#   counts_fn             Input counts table, txt format. Genes should be in raw counts.

#   design_fn             Input design table, must be a txt file. First column is condition, second column in contrast. Different contrasts are seperated by ","

# optional arguments:
#   -h, --help            show this help message and exit

#   --output [OUTPUT], -o [OUTPUT]
#     Output path. Default = "./results"

#   --gsea_ref [GSEA_REF], -r [GSEA_REF] 
#     GSEA library path. Must specify full path. Default = "/project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt".

#   --gsea_plot [GSEA_PLOT], -p [GSEA_PLOT] 
#     If GSEA plots will be made. By default no plots are made. Revert with "T"

#   --species [SPECIES], -s [SPECIES]
#     species from which the gene symbol is in. {human,mouse}

#   --fccutoff [FCCUTOFF], -f [FCCUTOFF]
#     Log fold change cutoff for volcano plot. Default 1

#   --vpcutoff [VPCUTOFF], -vp [VPCUTOFF]
#     Adjusted p-value cutoff for volcano plot. Default 0.05

#   --pcutoff [PCUTOFF], -hp [PCUTOFF]
#     Adjusted p-value cutoff for heatmap. Default 0.05

#   --partialheatmap, -ph
#     Toggle to only show relevant samples in DEG heatmap. Default FALSE

#   --nocolclustering, -nc
#     Toggle to turn off column clustering in DEG heatmap. Default TRUE

#   --ssgsea, -sg         
#     Toggle to run ssGSEA for all samples. Default FALSE

#   --wgcna, -w         
#     Toggle to run WGCNA for all samples. Default FALSE

# Try this:
# Rscript DE.r ./example_data/example_expression.txt ./example_data/example_group.txt -o ./results 
#   -r /project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt -g T
#   -f 1 -hp 0.05 -vp 0.001 -ph

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("DESeq2","fgsea",'fastmatch','pheatmap','GSVA', 'WGCNA', "HybridMTest")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(packages, rownames(installed.packages())))
}

packages <- c("ggplot2", "dplyr",'tibble','argparse', 'readr', 'gprofiler2')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

suppressMessages(library('DESeq2'))
suppressMessages(library('ggplot2'))
suppressMessages(library('ggrepel'))
suppressMessages(library('fgsea'))
suppressMessages(library('dplyr'))
suppressMessages(library('tibble'))
suppressMessages(library('pheatmap'))
suppressMessages(library('argparse'))
suppressMessages(library('GSVA'))
suppressMessages(library('gprofiler2'))

# debug testing for using in Rstudio
# setwd('/project/shared/xiao_wang/software/rnaseqDE/')
# design = './example_data/example_group.txt'
# cts = './example_data/example_expression.txt'
# output_path = './example_results/'
# gsea_ref = "/project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt"
# fc_cutoff <- 1
# pval_cutoff <- 0.05
# v_pval_cutoff <- 0.001
# species <- 'human'
# partialheatmap <- T
# clustercol <- T

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
  default='/project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt',
  help='GSEA library path. Default = "/project/shared/xiao_wang/software/rnaseqDE/example_data/h.all.v7.0.symbols.gmt".')

parser$add_argument(
  '--gsea_plot','-g',action = 'store_true', default = FALSE,
  help='If GSEA plots will be made. By default no plots are made. ')

parser$add_argument(
  '--gsea_heatmap','-gh',action = 'store_true', default = FALSE,
  help='If GSEA heatmap will be made. By default no plots are made. ')

parser$add_argument(
  '--species','-s',nargs='?', default='human',
  help='Species from which the gene symbol is in {human, mouse}. Default human')

parser$add_argument(
  '--fccutoff','-f',nargs='?', default=1, type='double',
  help='Log fold change cutoff for volcano plot. Default 1')

parser$add_argument(
  '--vpcutoff','-vp',nargs='?', default=0.05, type='double',
  help='Adjusted p-value cutoff for volcano plot. Default 0.05')

parser$add_argument(
  '--pcutoff','-hp',nargs='?', default=0.05, type='double',
  help='Adjusted p-value cutoff for heatmap. Default 0.05')
  
parser$add_argument(
  '--partialheatmap','-ph',action = 'store_true', default = FALSE,
  help = 'Toggle to only show relevant samples in DEG heatmap. Default FALSE')

parser$add_argument(
    '--nocolclustering', '-nc', action = 'store_false', default = TRUE,
    help = 'Toggle to turn off column clustering in DEG heatmap. Default TRUE')

parser$add_argument(
    '--ssgsea', '-sg', action = 'store_true', default = FALSE,
    help = 'Toggle to run ssGSEA for all samples. Default FALSE')

parser$add_argument(
    '--wgcna', '-w', action = 'store_true', default = FALSE,
    help = 'Toggle to run WGCNA for all samples. Default FALSE')

args <- parser$parse_args()
cts <- args$counts_fn
design <- args$design_fn
output_path <- args$output
gsea_ref <- args$gsea_ref
gsea_plot <- args$gsea_plot
gsea_heatmap <- args$gsea_heatmap
fc_cutoff <- args$fccutoff
pval_cutoff <- args$pcutoff
v_pval_cutoff <- args$vpcutoff
species <- args$species
partialheatmap <- args$partialheatmap
clustercol <- args$nocolclustering
ssgsea <- args$ssgsea
wgcna = args$wgcna

# --------
# Preprocessing input data

# gene names must not have repeats
design <- read.table(design,stringsAsFactors = F,header=T, sep='\t',row.names = 1)
design=design[order(design$Group),] # by wtwt5237
cts <- read.table(cts, stringsAsFactors = F,header=T, sep='\t', row.names = 1, check.names = F)

# Make output file and set path to it
dir.create(file.path(output_path), showWarnings = FALSE)
setwd(file.path(output_path))

# Setting up contrasts
contrast_groups <- unique(design)
analysis = list()
j=1
for (i in 1:dim(contrast_groups)[1]) {
  t <- contrast_groups[i,1]
  if (length(contrast_groups[i,2]) == 0) {next}
  refs <- strsplit(contrast_groups[i,2],',')[[1]]
  for (ref in refs) {
    analysis[[j]] <- list(t,ref)
    j <- j + 1
  }
} 

# Align dataset with design
design <- subset(design,row.names(design) %in% colnames(cts))
cts <- cts[,row.names(design)]

# Account for batch effects in design
if ("Batch" %in% colnames(design)){
  colnames(design) <- c('condition','Contrasts','Batch')
  dds <- DESeqDataSetFromMatrix(
    countData = cts, colData = design, design = ~ condition + Batch)
} else {
  colnames(design) <- c('condition','Contrasts')
  dds <- DESeqDataSetFromMatrix(
    countData = cts, colData = design, design = ~ condition)
}
# Simply checking for not expressed genes, filter at sum of counts of all
# samples >= 10
keep <- rowSums(counts(dds)) >= 2
dds = dds[keep,]

# PCA based on DESeq2 normalized counts
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData['Name'] = row.names(pcaData)
ggplot(
  pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(aes(label = Name), size = 4)
ggsave('PCA_plot.pdf',height = 5, width = 5)

# Distance plot
corr <- cor(assay(vsd), method='pearson')
pheatmap(
  corr,annotation_col = design['condition'],show_rownames=F,cluster_cols = T,
  filename = "Sample distance heatmap.pdf",sep = '')

# Write normalized count to results.
export_counts <- assay(vsd)
write.table(export_counts,'deseq2 vst counts.csv', sep=',')

# ssGSEA
if (ssgsea == T) {
  ssgsea_ref = gmtPathways(gsea_ref)
  ssgsea_res = gsva(export_counts,ssgsea_ref,method='ssgsea')
  pheatmap(
    ssgsea_res,annotation_col = design['condition'],cluster_cols = F, 
    height = 0.4*length(ssgsea_ref), width=0.4*dim(export_counts)[2],
    filename = "ssGSEA heatmap.pdf",sep = '', scale = 'row')
    write.table(ssgsea_res,'ssGSEA_results.csv', sep=',')
}

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
  data2voc = res[!is.na(res$padj),]
    data2voc$sig = ifelse(
        (abs(data2voc$log2FoldChange)>=fc_cutoff) & (data2voc$padj<=v_pval_cutoff),
        "Significant Up", "Non-significant")
    data2voc[
        (data2voc$sig == "Significant Up") & (data2voc$log2FoldChange<=-fc_cutoff),"sig"
        ] = "Significant Down"
  data2voc = as.data.frame(data2voc)
  data2voc$padj = -log10(data2voc$padj)
  cols=c("Significant Up"='red', "Significant Down"='blue',"Non-significant"='grey')
  # Capping logFC magnitude to 10
  data2voc$log2FoldChange = ifelse(data2voc$log2FoldChange > 10, 10, data2voc$log2FoldChange)
  data2voc$log2FoldChange = ifelse(data2voc$log2FoldChange < -10, -10, data2voc$log2FoldChange)
  
  data2voc = data2voc[order(abs(data2voc$log2FoldChange),decreasing = T),]
  keep = row.names(data2voc)[data2voc$sig %in% c('Significant Up','Significant Down')][1:20]
  keep = keep[!is.na(keep)]
  data2voc$Gene = NA
  data2voc[keep,"Gene"] = keep

  ggplot(
    data2voc,
    aes(log2FoldChange,padj,label=Gene))+scale_shape_discrete(solid=T)+geom_point(mapping=aes(col=sig),shape=16,size=1)+
    geom_line(aes(x=-fc_cutoff),linetype="dotted",color='gold')+
    geom_line(aes(x=fc_cutoff),linetype="dotted",color='gold')+
    geom_line(aes(y=-log10(v_pval_cutoff)),linetype="dotted",color='gold')+
    scale_color_manual(values=cols)+geom_text_repel(force=2,segment.size=0.25)+theme_bw()
  
  ggsave(paste(output_prefix,"/DEG_Volcano_plot.pdf",sep = '')) # by wtwt5237
  
  
  
  # Heatmap
  res_heatmap = res[res$padj<=pval_cutoff & !is.na(res$padj),] # by wtwt5237
  res_up = res_heatmap[res_heatmap$log2FoldChange > fc_cutoff, ]
  res_dn = res_heatmap[res_heatmap$log2FoldChange < -fc_cutoff, ]
  # by wtwt5237 -start
  tops = rownames(res_heatmap[order(abs(res_up$log2FoldChange), decreasing = T),])[1:min(100,dim(res_up)[1])]
  bots = rownames(res_heatmap[order(abs(res_dn$log2FoldChange), decreasing = T),])[1:min(100,dim(res_dn)[1])]
  top_genes = append(tops, bots)
  if (length(tops) <= 1) {
    next
  }
  if (partialheatmap == T) {
  heatmap_samples = row.names(design)[design$condition %in% c(target,ref)]
  } else {heatmap_samples = row.names(design)
  }
  heatmap_mats = export_counts[top_genes,heatmap_samples]
  # by wtwt5237 - end
  pheatmap(
    heatmap_mats,annotation_col = design['condition'],show_rownames=T, # by wtwt5237
    scale = 'row',cluster_cols = clustercol, height = 15,fontsize=6, # by wtwt5237 
    filename = paste(output_prefix,"/DEG_heatmap.pdf",sep = '')) # by wtwt5237
  
  # ========
  # GSEA analysis
  # pre-processing Deseq results
  res$gene = row.names(res)
  res2 = res[,c("gene","stat")]
  
  # If input is in mouse symbols, do id conversion before GSEA.
  if (species == 'mouse') {
    m2h = read.table(
      '/project/shared/xiao_wang/software/rnaseqDE/script/M2H_symbol_conversion.txt',
      stringsAsFactors = F,header=T)
    res2 = merge(as.data.frame(res2),m2h, by.x = 'gene', by.y = 'Symbol_x')
    res2 = aggregate(res2$stat, list(res2$Human_symbol), median)
    colnames(res2) = c('gene','stat')
  }
  
  ranks = deframe(res2)
  ranks = ranks[!is.na(ranks)]
  pathway_ref = gmtPathways(gsea_ref)
  
  # GSEA PreRank using LFC log2FoldChange in DESeq results, which is basically logFC.
  fgseaRes_kegg <- fgseaMultilevel(pathways=pathway_ref, stats=ranks)
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
          title="GSEA") +
    theme_minimal()
  ggsave(paste(output_prefix, "/GSEA.pdf",sep=""), width=9,height=16) # by wtwt5237
  
  # Make all GSEA plots
  if (gsea_plot == T){
    if(! dir.exists(paste(output_prefix,"GSEA_plots",sep=""))){
      dir.create(paste(output_prefix,"GSEA_plots",sep=""))
    }
    
    for (pathway_name in fgseaResTidy_c$pathway){
      plotEnrichment(pathway_ref[[pathway_name]], ranks)
      ggsave(paste(output_prefix,'/GSEA_plots/', pathway_name, '.pdf',sep='')) # by wtwt5237
    }
  }
  # Make all GSEA heatmaps
  if (gsea_heatmap == T){
    if(! dir.exists(paste(output_prefix,"GSEA_heatmaps",sep=""))){
      dir.create(paste(output_prefix,"GSEA_heatmaps",sep=""))
    }
    
  for (pathway_name in fgseaResTidy_c$pathway){
      pathway_genes = pathway_ref[pathway_name][[1]]
      pathway_genes = intersect(pathway_genes, row.names(export_counts))
      if (length(pathway_genes) >= 5) {
        heatmap_mat = export_counts[pathway_genes,heatmap_samples]
        heatmap_mat = heatmap_mat[apply(heatmap_mat,1,sd) > 0,]
        pheatmap(
          heatmap_mat,annotation_col = design['condition'],show_rownames=T, 
          scale = 'row',cluster_cols = clustercol, height = 15,fontsize=8, 
          filename = paste(output_prefix,"/GSEA_heatmaps/",pathway_name,".pdf",sep = '')
        )
      }
    }
  }
}

warnings() # by wtwt5237

# WGCNA
if (wgcna == T) {
  # 0. Setup
  suppressMessages(library('WGCNA'))
  datExpr = data.frame(t(export_counts))
  names(datExpr) = row.names(export_counts)
  dimnames(datExpr)[[1]] = names(data.frame(row.names(export_counts)))
  gsg = goodSamplesGenes(datExpr, verbose = 3)
  dir.create("WGCNA_Results")

  # 1. Scale Independence & Mean Connectivity
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  pdf(file = "WGCNA_Results/Scale Independence and Mean Connectivity.pdf")
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale Independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="blue");
  abline(h=0.90,col="blue")

  plot(sft$fitIndices[,1], sft$fitIndices[,5],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean Connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")
  dev.off()

  # 2. Cluster Dendrogram
  ysft <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  best_power = Position(function(x) x > 0.9, ysft)

  cor <- WGCNA::cor
  datExpr[] <- lapply(datExpr, as.numeric)
  net = blockwiseModules(datExpr, power = best_power,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0.01, mergeCutHeight = 0.1,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 3, nThreads = 16, deepSplit = 4)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  cor<-stats::cor

  sizeGrWindow(12, 9)
  mergedColors = labels2colors(net$colors)
  pdf(file = "WGCNA_Results/Cluster Dendrogram.pdf")
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  # 3. Average Expression Heatmap
  # reference <- table(net$colors)
  avgexp <- moduleEigengenes(datExpr, moduleColors)$averageExpr
  avgexp <- t(avgexp)
  avgexp = avgexp[,row.names((design))]

  suppressMessages(library('HybridMTest'))
  anova_res = row.oneway.anova(avgexp,design$Group)
  keep = anova_res$pval<=0.01
  avgexp = avgexp[keep,]

  pheatmap(
    avgexp, height = 16, fontsize = 8, cluster_cols = F,
    annotation_col = subset(design,select = 'Group'),
    scale='column',
    file = "WGCNA_Results/Average Expression Heatmap.pdf",
  )
  # 4. Output significant modules list
  sig_modules = row.names(avgexp)
  sig_modules = substring(sig_modules, first = 3, last = 99)
  geneInfo0 = data.frame(
    Genes = names(datExpr), 
    moduleLabel = moduleLabels, moduleColor = moduleColors)
  geneOrder = order(geneInfo0$moduleColor)
  geneInfo = geneInfo0[geneOrder, ]
  write.csv(geneInfo, "WGCNA_Results/geneInfo.csv")

  # 5. Significant module heatmap
  for (m in sig_modules) {
    genes = row.names(geneInfo[geneInfo$moduleColor == m,])
    heatmap_mats <- export_counts[genes, ]
    heatmap_mats <- na.omit(heatmap_mats)
    pheatmap(
      heatmap_mats, height = 16, fontsize = 5, 
      annotation_col = subset(design,select = 'Group'), cluster_cols = F,
      filename = paste("WGCNA_Results/Module_Heatmaps/", m, "_heatmap.pdf",sep = ''),
    )
  }

  # 6. Eigengene Heatmap & Dendrogram
  sig_modules <- paste0("ME", sig_modules)
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MET = orderMEs(MEs)
  cMET <- colnames(MET)
  cmatch <- match(sig_modules, cMET)
  sigMET <- MET[, cmatch]

  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  pdf(file = "WGCNA_Results/Eigengene Adjacency Heatmap.pdf")
  plotEigengeneNetworks(
    sigMET, "Eigengene Adjacency Heatmap", marDendro = c(0,4,1,2), 
    marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
  dev.off()

  # 7. Functional Enrichment for Significant Modules
  xli <- names(datExpr)
  genelist <- gconvert(xli, organism = "hsapiens", target = "ENTREZGENE_ACC", filter_na = T)
  genelist <- genelist[,!(names(genelist) %in% c("input_number", "target_number", "name", "description", "namespace"))]
  genelist <- genelist[- grep("nan", genelist$target),]
  entrezid <- genelist$target

  GOenr = GOenrichmentAnalysis(moduleColors, entrezid, organism = "human", nBestP = 10)
  tab = GOenr$bestPTerms[[4]]$enrichment

  sig_modules1 = sub("ME", "", sig_modules)
  sigtab <- tab[tab$module %in% sig_modules1, ]
  write.table(sigtab, file = "WGCNA_Results/GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
}