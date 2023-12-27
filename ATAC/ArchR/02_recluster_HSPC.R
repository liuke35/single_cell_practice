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

addArchRThreads(threads = 4)

subgroup <- "HSPC"
wd <- sprintf("~/Desktop/ATAC/202312/subclustered_%s", subgroup)
atac_sample_cmap = readRDS('~/Desktop/ATAC/202312/01_data_integration/samp_cmap.rds')
disease_cmap <- readRDS('~/Desktop/ATAC/202312/01_data_integration/disease_cmap.rds')
#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 1.0

###

atac_proj <- loadArchRProject(wd, force=TRUE)
library(zellkonverter)
rna_proj <- readH5AD('../HSPC_cell_raw.h5ad',reader = "python")

#atac_proj$group_info2 <- atac_proj$group_info
#atac_proj$group_info2[which(atac_proj$Sample == "ST_MA")] = 'aps_neg'
#atac_proj$group_info2[which(atac_proj$Sample == "ST_WR")] = 'aps_pos'

##########################################################################################
# Re-cluster subclustered ArchR project
##########################################################################################

atac_proj <- addIterativeLSI(
  ArchRProj = atac_proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  clusterParams = list(resolution=c(2), sampleCells=10000, n.start=10),
  sampleCellsPre = 30000,
  varFeatures = 25000,
  dimsToUse = 1:50,
  force = TRUE
)

# (Batch correct sample effects)
atac_proj <- addHarmony(
  ArchRProj = atac_proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  dimsToUse = 2:30,
  groupBy = "Sample",
  force = T
)

library(batchelor)
fMNN.res <- mnnCorrect(t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'HC_1'),1:50]),
                       t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'HC_2'),1:50]),
                       t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'ST_4'),1:50]),
                       t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'ST_5'),1:50]),
                       t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'ST_6'),1:50]),
                       t(atac_proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(atac_proj$Sample == 'ST_8'),1:50]))

atac_proj@reducedDims[['LSI_MNN']] <- SimpleList(
  matDR = t(fMNN.res@assays@data$corrected), 
  params = NA,
  date = Sys.time(),
  scaleDims = NA, #Do not scale dims after
  corToDepth = NA
)

# Identify Clusters from Iterative LSI
atac_proj <- addClusters(
  input = atac_proj,
  reducedDims = "LSI_MNN",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.3, 
  force = TRUE
)
#plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Clusters", embedding = "mnnUMAP")

set.seed(1)
atac_proj <- addUMAP(
  ArchRProj = atac_proj, 
  reducedDims = "LSI_MNN", 
  name = "UMAP", 
  nNeighbors = 35, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)

# Reorder clusters so they are sorted by cluster size
atac_proj <- relabelClusters(atac_proj)

# Make various cluster plots:
atac_proj <- addImputeWeights(atac_proj)
atac_proj <- visualizeClustering(atac_proj, pointSize=pointSize, embedding = 'UMAP',
                                 sampleCmap=atac_sample_cmap, diseaseCmap=disease_cmap)

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")


markerGenes <- c("CD34","AVP",
                 "DNTT","MME","VPREB1","RAG1",
                 "PAX5","JCHAIN",
                 "SPINK2","RUNX3","MKI67","CLC",
                 "GYPA","PF4","GATA1","GATA2","HBD",
                 "AZU1","ELANE","MPO")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(hm, name = "filtered-GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = atac_proj, addDOC = FALSE)

# Marker gene imputation with Magic:
p <- plotEmbedding(
  ArchRProj = atac_proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(atac_proj), 
  plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, width = 5, height = 5)

# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  tileSize=250,
  minCells=25
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, width = 5, height = 5)

####not
atac_proj <- addGeneIntegrationMatrix(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "LSI_MNN",
  seRNA = rna_proj, # Can be a seurat object
  addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
  force = TRUE,
  groupRNA = "cell_type", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
  nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
  nameGroup = "NamedClust_RNA", #Name of column where group from scRNA is matched to each cell
  nameScore = "predictedScore" #Name of column where prediction score from scRNA
)

pal <- paletteDiscrete(values = rna_proj$cell_type)

p_pred <- plotEmbedding(
  atac_proj, 
  colorBy = "cellColData", 
  name = "NamedClust_RNA", 
  embedding = "UMAP",
  pal = pal
)
###

# Relabel clusters:
plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "NamedClust", embedding = "UMAP")

nclust <- length(unique(atac_proj$Clusters))
fineClust <- c()
fineClust[1] <- 'Erythroid'
fineClust[2] <- 'Erythroid'
fineClust[3] <- 'MEP'
fineClust[4] <- 'NeuP'
fineClust[5] <- 'EBMP'
fineClust[6] <- 'HSC_MPP'
fineClust[7] <- 'Erythroid'
names(fineClust) <- atac_proj$Clusters %>% getFreqs() %>% names()
atac_proj$FineClust <- fineClust[atac_proj$Clusters] %>% unname

plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "FineClust", embedding = "UMAP")
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")
barwidth=0.9
namedClustBySamp <- fractionXbyY(atac_proj$FineClust, atac_proj$Sample, add_total=TRUE, xname="FineClust", yname="Sample")
pdf(paste0(plotDir, "/sampleByFineClustBarPlot.pdf"))
print(stackedBarPlot(namedClustBySamp, cmap=atac_sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

#atac_proj$group_celltype <- paste0(atac_proj$FineClust,'_',atac_proj$group_info)
saveArchRProject(atac_proj)

####differential analysis
#gene score matrix
library(ggrepel)

markers_df <- getMarkerFeatures(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "group_info",
  useGroups = "ST",
  bgdGroups = "HC",
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
  labs(title=" ST vs HC differential gene score") +
  geom_text_repel(aes(label = Label),size=2.5,color="black",direction="both",min.segment.length = 0.05,
                  segment.alpha=0.1,max.overlaps =30,nudge_x = 0.2,nudge_y=0.2) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black")  
  )

#do volcano plot for ST vs HC in each cell type to see differential TF motfs
library(BSgenome.Hsapiens.UCSC.hg38)

atac_proj <- addGroupCoverages(ArchRProj=atac_proj, groupBy="NamedClust")

pathToMacs2 <- '/Users/keliu/anaconda3/bin/macs2'

atac_proj <- addReproduciblePeakSet(
  ArchRProj = atac_proj, 
  groupBy = "FineClust", 
  pathToMacs2 = pathToMacs2
)

# Add Peak Matrix
atac_proj <- addPeakMatrix(ArchRProj = atac_proj)

atac_proj <- addMotifAnnotations(ArchRProj = atac_proj, motifSet = "JASPAR2020", annoName = "Motif")
atac_proj <- addBgdPeaks(atac_proj)
atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(atac_proj, name = "MotifMatrix", plot = TRUE)
saveArchRProject(atac_proj)

markersTF_df <- getMarkerFeatures(
  ArchRProj = atac_proj, 
  useMatrix = "MotifMatrix", 
  groupBy = "group_info",
  useGroups = "ST",
  bgdGroups = "HC",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

deTF <- as.data.frame(getMarkers(markersTF_df, cutOff = "FDR <= 1"))# 

deTF$abs_avg <- abs(deTF$MeanDiff)
#de$diffexpressed <- ifelse(de$Log2FC>0,'up','down')
deTF$diffexpressed <- "NO"
deTF$diffexpressed[deTF$MeanDiff > 0.01 & deTF$FDR < 0.05] <- "up"
deTF$diffexpressed[deTF$MeanDiff < -0.01 & deTF$FDR < 0.05] <- "down"
#de$TF <- str_split_fixed(de$name,'_',2)[,1]
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
  labs(title=" ST vs HC differential TF score") +
  geom_text_repel(aes(label = Label),size=2.5,color="black",direction="both",min.segment.length = 0.05,
                  segment.alpha=0.1,max.overlaps =30,nudge_x = 0.01,nudge_y=0.01) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black")  
  )
