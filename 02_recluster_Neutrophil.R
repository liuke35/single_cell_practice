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

subgroup <- "Neutrophil"
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
rna_proj <- readH5AD('../Neu_cell_raw.h5ad',reader = "python")

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
  resolution = 0.5, 
  force = TRUE
)
plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

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


markerGenes <- c("CD34", "HLF", "AVP", "SPINK2","AZU1","ELANE",  # HSPC
                 "CLC","GZMB",
                 "CEBPE","ITGAM", "MMP9", "PGLYRP1", "ORM1", "CEACAM8","PADI4")

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
####

# Relabel clusters:
plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "NamedClust", embedding = "UMAP")

nclust <- length(unique(atac_proj$Clusters))
fineClust <- c()
fineClust[1] <- 'mature_Neutrophil'
fineClust[2] <- 'immature_Neutrophil'
fineClust[3] <- 'mature_Neutrophil'
fineClust[4] <- 'immature_Neutrophil'
fineClust[5] <- 'early_immature_Neutrophil'
fineClust[6] <- 'immature_Neutrophil'
fineClust[7] <- 'early_immature_Neutrophil'
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
