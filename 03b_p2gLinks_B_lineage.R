#!/usr/bin/env Rscript

#############################
# Processing scalp data
#############################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# Get additional functions, etc.:
source('~/Desktop/ATAC/analysis_function/archr_helper.R')
source('~/Desktop/ATAC/analysis_function/plotting_config.R')
source('~/Desktop/ATAC/analysis_function/misc_helper.R')
source('~/Desktop/ATAC/analysis_function/matrix_helper.R')
source('~/Desktop/ATAC/analysis_function/homemade_archr.R')
source('~/Desktop/ATAC/analysis_function/GO_wrappers.R')
# Set Threads to be used
addArchRThreads(threads = 4)

# set working directory
subgroup <- "B_lineage"
wd <- sprintf("~/Desktop/ATAC/202312/subclustered_%s", subgroup)
full_dir <- "~/Desktop/ATAC/202312/fine_clustered"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.0

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
#rna_proj <- readRDS(sprintf("/oak/stanford/groups/wjg/boberrey/hairATAC/results/scRNA_preprocessing/harmonized_subclustering/%s/%s.rds", subgroup, subgroup))

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# Colormaps
sample_cmap = readRDS('~/Desktop/ATAC/202312/01_data_integration/samp_cmap.rds')
disease_cmap <- readRDS('~/Desktop/ATAC/202312/01_data_integration/disease_cmap.rds')

###########################################################################################
# Do not proceed prior to calling peaks
###########################################################################################
library(BSgenome.Hsapiens.UCSC.hg38)

# Compute group coverages
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)

# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(atac_proj$FineClust), peakGR=full_peaks)
rm(full_proj); gc()

# Now add these peaks to the subproject and generate peak matrix
atac_proj <- addPeakSet(atac_proj, peakSet=peaks, force=TRUE)
atac_proj <- addPeakMatrix(atac_proj, force=TRUE)
atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

# Calculate coaccessibility
atac_proj <- addCoAccessibility(
  ArchRProj = atac_proj,
  reducedDims = "LSI_MNN"
)

# Calculate peak-to-gene links
atac_proj <- addPeak2GeneLinks(
  ArchRProj = atac_proj,
  reducedDims = "LSI_MNN"
)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

# Add deviations matrix
atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Save intermediate output
saveArchRProject(atac_proj)

p_p2g <- plotPeak2GeneHeatmap(ArchRProj = atac_proj, groupBy = "FineClust")
###########################################################################################

library(ggrepel)

atac_proj$group_celltype <- paste0(atac_proj$FineClust,'_',atac_proj$group_info)

for (i in names(table(atac_proj$FineClust))){
  markers_df <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "group_celltype",
    useGroups = paste0(i,"_ST"),
    bgdGroups = paste0(i,"_HC"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  de <- as.data.frame(getMarkers(markers_df, cutOff = "FDR <= 1"))# 
  
  de$abs_avg <- abs(de$Log2FC)
  #de$diffexpressed <- ifelse(de$Log2FC>0,'up','down')
  de$diffexpressed <- "NO"
  de$diffexpressed[de$Log2FC > 0.5 & de$FDR < 0.05] <- "up"
  de$diffexpressed[de$Log2FC < -0.5 & de$FDR < 0.05] <- "down"
  #de$TF <- str_split_fixed(de$name,'_',2)[,1]
  de <- de[order(de$FDR,-de$abs_avg),]
  up.genes <- head(de$name[which(de$diffexpressed == "up")], 10)
  
  down.genes <- head(de$name[which(de$diffexpressed == "down")], 10)
  
  de$Label = ""
  deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
  de$Label[match(deg.top10.genes, de$name)] <- deg.top10.genes
  p <- ggplot(data=de, aes(x=Log2FC, y=-log10(FDR), col = diffexpressed,label=Label)) +
    scale_color_manual(values=c("blue", "grey", "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="black",linetype="dashed",alpha = 0.5) +
    geom_hline(yintercept=-log10(0.05), col="black",linetype="dashed",alpha = 0.5)+
    geom_point(size = .1) + 
    theme_minimal() +
    labs(title=paste0("ST vs HC differential gene score in ",i)) +
    geom_text_repel(aes(label = Label),size=2.5,color="black",direction="both",min.segment.length = 0.05,
                    segment.alpha=0.1,max.overlaps =30,nudge_x = 0.2,nudge_y=0.2) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.line = element_line(color = "black")  
    )
  ggsave(paste0(plotDir,"/DEG_score_ST_ ",i,".pdf"),p,width = 10,height = 10)
}

#do volcano plot for ST vs HC in each cell type to see differential TF motfs
plotVarDev <- getVarDeviations(atac_proj, name = "MotifMatrix", plot = TRUE)

for (i in names(table(atac_proj$FineClust))){
  markersTF_df <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "MotifMatrix", 
    groupBy = "group_celltype",
    useGroups = paste0(i,"_ST"),
    bgdGroups = paste0(i,"_HC"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  deTF <- as.data.frame(getMarkers(markersTF_df, cutOff = "FDR <= 1"))# 
  
  deTF$abs_avg <- abs(deTF$MeanDiff)
  deTF$diffexpressed <- "NO"
  deTF$diffexpressed[deTF$MeanDiff > 0.01 & deTF$FDR < 0.05] <- "up"
  deTF$diffexpressed[deTF$MeanDiff < -0.01 & deTF$FDR < 0.05] <- "down"
  deTF <- deTF[order(deTF$FDR,-deTF$abs_avg),]
  up.genes <- head(deTF$name[which(deTF$diffexpressed == "up")], 10)
  
  down.genes <- head(deTF$name[which(deTF$diffexpressed == "down")], 10)
  
  deTF$Label = ""
  deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
  deTF$Label[match(deg.top10.genes, deTF$name)] <- deg.top10.genes
  p <- ggplot(data=deTF, aes(x=MeanDiff, y=-log10(FDR), col = diffexpressed,label=Label)) +
    scale_color_manual(values=c("blue", "grey", "red")) +
    geom_vline(xintercept=c(-0.01, 0.01), col="black",linetype="dashed",alpha = 0.5) +
    geom_hline(yintercept=-log10(0.05), col="black",linetype="dashed",alpha = 0.5)+
    geom_point(size = 1) + 
    theme_minimal() +
    labs(title=paste0("ST vs HC differential TF score in ",i)) +
    geom_text_repel(aes(label = Label),size=2.5,color="black",direction="both",min.segment.length = 0.05,
                    segment.alpha=0.1,max.overlaps =30,nudge_x = 0.01,nudge_y=0.01) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.line = element_line(color = "black")  
    )
  ggsave(paste0(plotDir,"/DTF_score_ST_ ",i,".pdf"),p,width = 10,height = 10)
}


##########################################################################################
# Identify Correlated TF Motifs and TF Gene Score/Expression
##########################################################################################

seGroupMotif <- getGroupSE(ArchRProj = atac_proj, useMatrix = "MotifMatrix", groupBy = "FineClust")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = atac_proj,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_MNN"
)

corGIM_MM <- correlateMatrices(
  ArchRProj = atac_proj,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_MNN"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )+ ggrepel::geom_text_repel(
    data=data.frame(corGSM_MM)[data.frame(corGSM_MM)$TFRegulator=="YES",], aes(x=cor, y=maxDelta, label=GeneScoreMatrix_name), 
    size=2,
    point.padding=0, # additional pading around each point
    box.padding=0.5,
    min.segment.length=0, # draw all line segments
    max.overlaps=Inf, # draw all labels
    #nudge_x = 2,
    color="black"
  ) 

pdf(paste0(plotDir, "/corGSM_MM_posTFregulators.pdf"), width=10, height=10)
p
dev.off()

corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )+ ggrepel::geom_text_repel(
    data=data.frame(corGIM_MM)[data.frame(corGIM_MM)$TFRegulator=="YES",], aes(x=cor, y=maxDelta, label=GeneIntegrationMatrix_name), 
    size=2,
    point.padding=0, # additional pading around each point
    box.padding=0.5,
    min.segment.length=0, # draw all line segments
    max.overlaps=Inf, # draw all labels
    #nudge_x = 2,
    color="black"
  ) 

pdf(paste0(plotDir, "/corGIM_MM_posTFregulators.pdf"), width=10, height=10)
p
dev.off()

##########################################################################################
# Identify regulatory targets of TFs 
##########################################################################################
motifMatrix <- getMatrixFromProject(atac_proj, useMatrix="MotifMatrix")
deviationsMatrix <- assays(motifMatrix)$deviations

# GeneIntegration Matrix: (rows gene names x cols cell names)
GIMatrix <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
GImat <- assays(GIMatrix)$GeneIntegrationMatrix
rownames(GImat) <- rowData(GIMatrix)$name
GImat <- as(GImat[Matrix::rowSums(GImat) > 0,], "sparseMatrix") # Remove unexpressed genes

# Use only motifs that are 'TFRegulators' as determined by analysis above
GSMreg <- rownames(motifMatrix)[corGSM_MM[corGSM_MM$TFRegulator == "YES",]$MotifMatrix_idx]
GIMreg <- rownames(motifMatrix)[corGIM_MM[corGIM_MM$TFRegulator == "YES",]$MotifMatrix_idx]
regulators <- unique(c(GSMreg, GIMreg))
regulators <- c(
  "CTCF_177", "CTCFL_198", "ZBTB3_269", "HOXD8_537", "AR_689", "FLI1_337", 
  "PITX2_504", "ETV1_320", "ETV5_347", "ETV4_345",
  regulators
)
deviationsMatrix <- deviationsMatrix[regulators,]

# Identify pseudobulks for performing matrix correlations
knn_groups <- getLowOverlapAggregates(atac_proj, target.agg=500, k=100, 
                                      overlapCutoff=0.8, dimReduc="LSI_MNN")

kgrps <- unique(knn_groups$group)

# GeneIntegrationMatrix
GIMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(GImat[,use_cells])
}) %>% do.call(cbind,.)
colnames(GIMatPsB) <- kgrps

# In rare instances, we can get pseudo-bulked genes that have zero averages
GIMatPsB <- GIMatPsB[Matrix::rowSums(GIMatPsB) > 0,]


# DeviationsMatrix
DevMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(deviationsMatrix[,use_cells])
}) %>% do.call(cbind,.)
colnames(DevMatPsB) <- kgrps

# Perform chromVAR deviations to Integrated RNA correlation analysis:
start <- Sys.time()
geneCorMat <- cor2Matrices(DevMatPsB, GIMatPsB)
colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
end <- Sys.time()
message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))

allGenes <- rownames(GIMatPsB) %>% sort() # Already filtered to only expressed genes

# Get locations of motifs of interest:
motifPositions <- getPositions(atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# Get peak to gene GR
corrCutoff=0.45
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)

calculateLinkageScore <- function(motifLocs, p2gGR){
  # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
  ###################################
  # For a given gene, the LS = sum(corr peak R2 * motifScore)
  ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
  olGenes <- p2gGR[to(ol)]
  olGenes$motifScore <- motifLocs[from(ol)]$score
  olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
  LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
  LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
  LSdf$rank <- 1:nrow(LSdf)
  return(LSdf)
}

calculateMotifEnrichment <- function(motifLocs, p2gGR){
  # Calculate Motif enrichment per gene
  ###################################
  # For a given gene, calculate the hypergeometric enrichment of motifs in 
  # linked peaks (generally will be underpowered)
  motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
  m <- length(motifP2G) # Number of possible successes in background
  n <- length(p2gGR) - m # Number of non-successes in background
  
  motifLinks <- motifP2G$symbol %>% getFreqs()
  allLinks <- p2gGR$symbol %>% getFreqs()
  df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
  df$motifLinks[is.na(df$motifLinks)] <- 0
  df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
  df <- df[order(df$mLog10pval, decreasing=TRUE),]
  df$symbol <- rownames(df)
  return(df)
}


# plot all TF regulators
regPlotDir <- paste0(plotDir, "/TFregulatorPlots")
dir.create(regPlotDir, showWarnings = FALSE, recursive = TRUE)

markerGenes <- c("CD34","AVP",
                 "DNTT","MME","VPREB1","RAG1",
                 "PAX5","MS4A1","JCHAIN","CD19",
                 "SPINK2","RUNX3","IRF7","MKI67") %>% unique()

# Get list of genes we want to highlight (e.g. genes involved in HF development)
library(org.Hs.eg.db)
library(GO.db)

search_term <- "Interferon"  # 用你想搜索的字符串替换
matching_terms <- Term(GOTERM)  # 获取所有GO词条
matching_terms <- matching_terms[grep(search_term, matching_terms, ignore.case = TRUE)]  # 使用正则表达式搜索

go_id = GOID(GOTERM[Term(GOTERM) == "B cell activation"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
b_activation_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

go_id = GOID(GOTERM[Term(GOTERM) == "response to type I interferon"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
interferon1_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

go_id = GOID(GOTERM[Term(GOTERM) == "response to type II interferon"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
interferon2_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

go_id = GOID(GOTERM[Term(GOTERM) == "response to type III interferon"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)
interferon3_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()

markerGenes <- c(b_activation_genes, interferon1_genes, interferon2_genes, interferon3_genes, markerGenes) %>% unique() %>% sort()

# Store results for each TF
res_list <- list()

for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- geneCorMat[geneCorMat$motifName == motif,]
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.8)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes),c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0(regPlotDir, sprintf("/regulatory_targets_%s.rds", motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
    #+ geom_point(size = 2)
    + ggrastr::geom_point_rast(size=2)
    + ggrepel::geom_text_repel(
      data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
      #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
      size=2,
      point.padding=0, # additional pading around each point
      box.padding=0.5,
      min.segment.length=0, # draw all line segments
      max.overlaps=Inf, # draw all labels
      #nudge_x = 2,
      color="black"
    ) 
    + geom_vline(xintercept=0, lty="dashed") 
    + geom_vline(xintercept=corr_window, lty="dashed", color="red")
    + geom_vline(xintercept=-corr_window, lty="dashed", color="red")
    + geom_hline(yintercept=LS_window, lty="dashed", color="red")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1.0,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle=90, hjust=1))
    + ylab("Linkage Score") 
    + xlab("Motif Correlation to Gene") 
    + scale_color_gradientn(colors=cmaps_BOR$zissou, limits=c(0, maxPval))
    + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
    + scale_x_continuous(limits = c(-0.85, 0.955)) # Force plot limits
    + ggtitle(sprintf("%s putative targets", motif_short))
  )
  # Positively regulated genes:
  upGO <- rbind(
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="MF")
  )
  upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
                             barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(pos_top_genes)), enrichLimits=c(0, 6))
  # Negatively regulated genes:
  downGO <- rbind(
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="MF")
  )
  downGO <- downGO[order(as.numeric(downGO$pvalue), decreasing=FALSE),]
  down_go_plot <- topGObarPlot(downGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
                               barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(neg_top_genes)), enrichLimits=c(0, 6))
  pdf(paste0(regPlotDir, sprintf("/%s_LS.pdf", motif_short)), width=8, height=6)
  print(p)
  print(up_go_plot)
  print(down_go_plot)
  dev.off()
}

###########################################################################################
# motif footprint 
###########################################################################################
fpdir <- "../B_fp_subset"

fp_proj <- subsetArchRProject(
  ArchRProj = atac_proj,
  cells = getCellNames(atac_proj),
  outputDirectory = fpdir,
  dropCells = TRUE,
  force = TRUE
)
saveArchRProject(fp_proj)
rm(atac_proj);gc()
fp_proj$group_celltype <- paste0(fp_proj$FineClust,'_',fp_proj$group_info)

# Delete stuff we don't need to duplicate
unlink(paste0(fpdir, "/Plots/*"), recursive = TRUE)
unlink(paste0(fpdir, "/GroupCoverages"), recursive = TRUE)
unlink(paste0(fpdir, "/PeakCalls"), recursive = TRUE)
unlink(paste0(fpdir, "/*_filtered_barcodes.txt"), recursive = TRUE)

motifPositions <- getPositions(fp_proj)

motifs <- c("CTCF", "CTCFL", "ETV","ELK", "ZBTB3", "HOXD8", "PITX2","FLI1")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "ZBTB33_255"] %>% unique()

fp_proj <- addGroupCoverages(
  ArchRProj=fp_proj, 
  groupBy="group_info", 
  minCells = 50, 
  force=TRUE
)

seFoot <- getFootprints(
  ArchRProj = fp_proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group_info"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = fp_proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-group",
  addDOC = FALSE,
  smoothWindow = 5
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = fp_proj, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias-group",
  addDOC = FALSE,
  smoothWindow = 5
)

fp_proj <- addGroupCoverages(
  ArchRProj=fp_proj, 
  groupBy="group_celltype", 
  minCells = 50, 
  force=TRUE
)

seFoot <- getFootprints(
  ArchRProj = fp_proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group_celltype"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = fp_proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-group_celltype",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = fp_proj, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias-group_celltype",
  addDOC = FALSE,
  smoothWindow = 5
)


