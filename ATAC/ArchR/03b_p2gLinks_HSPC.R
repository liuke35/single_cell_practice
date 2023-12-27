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
subgroup <- "HSPC"
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


###########################################################################################

library(ggrepel)

atac_proj$group_celltype <- paste0(atac_proj$FineClust,'_',atac_proj$group_info)

p_p2g <- plotPeak2GeneHeatmap(ArchRProj = atac_proj, groupBy = "FineClust")



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

###########################################################################################
# DAR
###########################################################################################
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

st_filelist <- list()
group_dap <- list()
dapDir <- paste0(plotDir, "/DAP")
dir.create(dapDir, showWarnings = FALSE, recursive = TRUE)

j=1
for (i in names(table(atac_proj$FineClust))){
  markersPeak_df <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "group_celltype",
    useGroups = paste0(i,"_ST"),
    bgdGroups = paste0(i,"_HC"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  deP <- as.data.frame(getMarkers(markersPeak_df, cutOff = "Log2FC >= 0.25"))
  group_dap[[i]] <-deP
  df <- deP[,c('seqnames','start','end')]
  colnames(df) <- c('chrom','start','end')
  st_filelist[[j]] <- paste0(dapDir,'/',i,'_','st_dap.bed')
  write.table(df,file = paste0(dapDir,'/',i,'_','st_dap.bed'),quote = F,sep = '\t', row.names = F, col.names = F)
  j = j+1
}

for (i in names(group_dap)){
  group_dap[[i]]$gene <- paste0(group_dap[[i]]$seqnames,'-',group_dap[[i]]$start,'-',group_dap[[i]]$end)
}

bed_list <- lapply(st_filelist,readPeakFile)
peakAnno_list <- lapply(bed_list,annotatePeak,tssRegion = c(-3000,3000),
                        TxDb = txdb,
                        addFlankGeneInfo = T, flankDistance = 5000,
                        annoDb = "org.Hs.eg.db")
names(peakAnno_list) <- names(group_dap)
peakAnno_df_list <- lapply(peakAnno_list,as.data.frame)
for (i in 1:length(peakAnno_df_list)){
  peakAnno_df_list[[i]]$chr_pos <- paste0(peakAnno_df_list[[i]]$seqnames,'-',
                                          peakAnno_df_list[[i]]$start-1,'-',
                                          peakAnno_df_list[[i]]$end)
  peakAnno_df_list[[i]] <- merge(peakAnno_df_list[[i]],
                                 group_dap[[i]],
                                 by.x = 'chr_pos', by.y = 'gene')
}

saveRDS(peakAnno_list,'st_peakAnno_list.rds')

p1 <- plotAnnoBar(peakAnno_list)
p2 <- plotDistToTSS(peakAnno_list)
ggsave(paste0(dapDir,'/','HSPC_DAP_st_annotation.pdf'),p1+p2,width = 10, height = 5)

genes <- lapply(peakAnno_list, function(i) 
  as.data.frame(i)$geneId)
vennplot(genes[2:4], by='Vennerable')

genes_pro <- lapply(peakAnno_df_list, function(i) 
  filter(i,annotation == "Promoter (<=1kb)"))
genes_pro <- lapply(genes_pro, function(i) 
  i$geneId)

cc = compareCluster(geneCluster   = genes_pro,
                    fun           = "enrichPathway",
                    pvalueCutoff  = 0.05)

dotplot(cc, showCategory=5)
####
for (i in names(table(atac_proj$FineClust))){
  markerTest <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "group_celltype",
    useGroups = paste0(i,"_ST"),
    bgdGroups = paste0(i,"_HC"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
  motif_plot(atac_proj = atac_proj,markerTest = markerTest,prefix = paste0("ST-vs-HC_",i))
}

###########################################################################################
# motif footprint 
###########################################################################################
fpdir <- "../HSPC_fp_subset"

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

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "RUNX1", "SPI1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]

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

plotGroups(
  ArchRProj = fp_proj,
  groupBy = "group_info",
  colorBy = "MotifMatrix",
  name = "z:RUNX1_733",
  imputeWeights = getImputeWeights(fp_proj),
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotGroups(
  ArchRProj = fp_proj,
  groupBy = "group_celltype",
  colorBy = "MotifMatrix",
  name = "z:RUNX1_733",
  imputeWeights = getImputeWeights(fp_proj),
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE
)

p <- plotEmbedding(
  ArchRProj = fp_proj, 
  colorBy = "GeneScoreMatrix", 
  name = c("GFI1","NFIB","MECOM","HOXB5"), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(fp_proj), 
  plotAs="points", size = 1.0
)
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = fp_proj, 
        addDOC = FALSE, width = 5, height = 5)






