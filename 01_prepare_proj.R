suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(ggrastr)
})

source('~/Desktop/ATAC/analysis_function/archr_helper.R')
source('~/Desktop/ATAC/analysis_function/plotting_config.R')
source('~/Desktop/ATAC/analysis_function/misc_helper.R')
source('~/Desktop/ATAC/analysis_function/matrix_helper.R')

addArchRThreads(threads = 6)

wd <- "~/Desktop/ATAC/202312/01_data_integration"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

inputFiles <- c("~/Desktop/ATAC/LZY-BM-ATAC/fragments.tsv.gz",
                "~/Desktop/ATAC/DX-BM-ATAC/fragments.tsv.gz",
                "~/Desktop/ATAC/WR-BM-ATAC/fragments.tsv.gz",
                "~/Desktop/ATAC/MA-BM-ATAC/fragments.tsv.gz",
                "~/Desktop/ATAC/LM-BM-ATAC/fragments.tsv.gz",
                "~/Desktop/ATAC/ZWY-BM-ATAC/fragments.tsv.gz")
names(inputFiles) <- c("HC_1","HC_2","ST_4","ST_5","ST_6","ST_8")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnno = geneAnno,
  genomeAnno = genomeAnno,
  minTSS = 0, # Don't filter at this point
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)

#Create ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "unfiltered_output"
)

# Remove Arrow files after copying
unlink(paste0(wd, "/*.arrow"))

# Now, identify likely cells:
identifyCells <- function(df, TSS_cutoff=6, nFrags_cutoff=2000, minTSS=5, minFrags=1000, maxG=4){
    # Identify likely cells based on gaussian mixture modelling.
    # Assumes that cells, chromatin debris, and other contaminants are derived from
    # distinct gaussians in the TSS x log10 nFrags space. Fit a mixture model to each sample
    # and retain only cells that are derived from a population with mean TSS and nFrags passing
    # cutoffs
    ####################################################################
    # df = data.frame of a single sample with columns of log10nFrags and TSSEnrichment
    # TSS_cutoff = the TSS cutoff that the mean of a generating gaussian must exceed
    # nFrags_cutoff = the log10nFrags cutoff that the mean of a generating gaussian must exceed
    # minTSS = a hard cutoff of minimum TSS for keeping cells, regardless of their generating gaussian
    # maxG = maximum number of generating gaussians allowed

    cellLabel <- "cell"
    notCellLabel <- "not_cell"

    if(nFrags_cutoff > 100){
        nFrags_cutoff <- log10(nFrags_cutoff)
        minFrags <- log10(minFrags)
    } 
    
    # Fit model
    set.seed(1)
    mod <- Mclust(df, G=2:maxG, modelNames="VVV")

    # Identify classifications that are likely cells
    means <- mod$parameters$mean

    # Identify the gaussian with the maximum TSS cutoff
    idents <- rep(notCellLabel, ncol(means))
    idents[which.max(means["TSSEnrichment",])] <- cellLabel

    names(idents) <- 1:ncol(means)

    # Now return classifications and uncertainties
    df$classification <- idents[mod$classification]
    df$classification[df$TSSEnrichment < minTSS] <- notCellLabel
    df$classification[df$nFrags < minFrags] <- notCellLabel
    df$cell_uncertainty <- NA
    df$cell_uncertainty[df$classification == cellLabel] <- mod$uncertainty[df$classification == cellLabel]
    return(list(results=df, model=mod))
}

# Run classification on all samples
minTSS <- 3
samples <- unique(proj$Sample)
cellData <- getCellColData(proj)
cellResults <- lapply(samples, function(x){
  df <- cellData[cellData$Sample == x,c("nFrags","TSSEnrichment")]
  df$log10nFrags <- log10(df$nFrags)
  df <- df[,c("log10nFrags","TSSEnrichment")]
  identifyCells(df, minTSS=minTSS)
  })
names(cellResults) <- samples

# Save models for future reference
saveRDS(cellResults, file = paste0(proj@projectMetadata$outputDirectory, "/cellFiltering.rds"))

# Plot filtering results
for(samp in samples){
    df <- as.data.frame(cellResults[[samp]]$results)
    cell_df <- df[df$classification == "cell",]
    non_cell_df <- df[df$classification != "cell",]

    xlims <- c(log10(500), log10(100000))
    ylims <- c(0, 18)
    # QC Fragments by TSS plot w/ filtered cells removed:
    p <- ggPoint(
        x = cell_df[,1], 
        y = cell_df[,2], 
        size = 1.5,
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = xlims,
        ylim = ylims,
        title = sprintf("%s droplets plotted", nrow(cell_df)),
        rastr = TRUE
    )
    # Add grey dots for non-cells
    p <- p + geom_point_rast(data=non_cell_df, aes(x=log10nFrags, y=TSSEnrichment), color="light grey", size=0.5)
    p <- p + geom_hline(yintercept = minTSS, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
    plotPDF(p, name = paste0(samp,"_EM_model_filtered_cells_TSS-vs-Frags.pdf"), ArchRProj = proj, addDOC = FALSE)
}

# Now, filter ATAC project to contain only cells
finalCellCalls <- lapply(cellResults, function(x) x$results) %>% do.call(rbind, .)
proj <- addCellColData(proj, data=finalCellCalls$classification, name="cellCall", cells=rownames(finalCellCalls), force=TRUE)
proj <- addCellColData(proj, data=finalCellCalls$cell_uncertainty, name="cellCallUncertainty", cells=rownames(finalCellCalls), force=TRUE)
group_info <-substr(proj$Sample,1,2)
proj$group_info <- group_info

#violin plot of unfiltered project
#HC_1  HC_2  ST_4  ST_5  ST_6  ST_8 
#13133 23470 16481 14593 63600 20442 
p1 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p2 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p3 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "BlacklistRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


##########################################################################################
# two plans of QC, the first is a harder filter from a 2023 Nature Genetics paper
# the second is a simpler filter that tutorial usually uses
##########################################################################################

# first method:Real cells pass QC filter 
realCells <- getCellNames(proj)[(proj$cellCall == "cell")]
subProj <- subsetArchRProject(proj, cells=realCells, 
    outputDirectory="filtered_output", dropCells=TRUE, force=TRUE)

# second method:Real cells pass QC filter 
idxPass <- which(proj$TSSEnrichment >= 10)
idxPass <- intersect(idxPass,which(proj$BlacklistRatio < 0.05))
idxPass <- intersect(idxPass,which(proj$nFrags > 3000))
idxPass <- intersect(idxPass,which(proj$NucleosomeRatio < 4))
naive_filter <- getCellNames(proj)[idxPass]
subProj_naive <- subsetArchRProject(proj, cells=naive_filter, 
                              outputDirectory="naive_filtered_output", dropCells=TRUE, force=TRUE)

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

subProj_naive <- addTileMatrix(subProj_naive, force=TRUE)
subProj_naive <- addGeneScoreMatrix(subProj_naive, force=TRUE)

# Add Infered Doublet Scores to ArchR project (~5-10 minutes)
subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)
subProj_naive <- addDoubletScores(subProj_naive, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)

########################################################################################
#in order to simplify codes, we assign subProj/subProj_naive to subProj to continue,
#remember to save the result of each method
########################################################################################
subProj <- subProj_naive
rm(subProj_naive)
###########

plotList <- list()
plotList[[1]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment"
)
plotList[[2]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "DoubletEnrichment"
)
plotPDF(plotList = plotList, name = "TSS-Doublet-Enrichment", width = 4, height = 4,  ArchRProj = subProj, addDOC = FALSE)

subProj <- filterDoublets(subProj, filterRatio = 1)
saveArchRProject(subProj)

##########################################################################################
# Reduced Dimensions and Clustering
##########################################################################################

subProj <- addIterativeLSI(
    ArchRProj = subProj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    sampleCellsPre = 15000,
    varFeatures = 50000, 
    dimsToUse = 1:50,
    force = TRUE
)

subProj <- addHarmony(
  ArchRProj = subProj,
  reducedDims = "IterativeLSI",
  dimsToUse = 2:30,
  name = "Harmony",
  groupBy = "Sample",
  force = T
)

library(batchelor)
fMNN.res <- mnnCorrect(t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'HC_1'),1:50]),
                       t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'HC_2'),1:50]),
                       t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'ST_4'),1:50]),
                       t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'ST_5'),1:50]),
                       t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'ST_6'),1:50]),
                       t(subProj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(subProj$Sample == 'ST_8'),1:50]))


subProj@reducedDims[['LSI_MNN']] <- SimpleList(
  matDR = t(fMNN.res@assays@data$corrected), 
  params = NA,
  date = Sys.time(),
  scaleDims = NA, #Do not scale dims after
  corToDepth = NA
)
# Identify Clusters from Iterative LSI
subProj <- addClusters(
    input = subProj,
    reducedDims = "LSI_MNN",
    dimsToUse = 2:30,
    method = "Seurat",
    name = "Clusters",
    maxClusters = 30,
    resolution = 0.6,
    force = TRUE
)

##########################################################################################
# Visualize Data
##########################################################################################

set.seed(1)
subProj <- addUMAP(
    ArchRProj = subProj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    dimsToUse = 2:30,
    nNeighbors = 50, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

subProj <- addUMAP(
  ArchRProj = subProj, 
  reducedDims = "Harmony", 
  name = "HarmonyUMAP", 
  dimsToUse = 2:30,
  nNeighbors = 50, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)

subProj <- addUMAP(
  ArchRProj = subProj, 
  reducedDims = "LSI_MNN", 
  name = "mnnUMAP", 
  dimsToUse = 2:30,
  nNeighbors = 50, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)
# Relabel clusters so they are sorted by cluster size
subProj <- relabelClusters(subProj)

subProj <- addImputeWeights(subProj)
samp_cmap <- c("#FF7D7D","#ba79b1","#896d47","#86908a","#e18a3a","#a45e44")
names(samp_cmap) <- names(table(subProj$Sample))
disease_cmap <- c("#d3a237","#4f794a")
names(disease_cmap) <- names(table(subProj$group_info))
saveRDS(samp_cmap,'samp_cmap.rds')
saveRDS(disease_cmap,'disease_cmap.rds')

# Make various cluster plots:
pointSize = 0.5
subProj <- visualizeClustering(subProj, pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)

p1 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Sample", embedding = "HarmonyUMAP")
p2 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Clusters", embedding = "HarmonyUMAP")

p3 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Sample", embedding = "mnnUMAP")
p4 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Clusters", embedding = "mnnUMAP")

plotPDF(p1,p2,p3,p4, name = "Plot-BatchCorrection-UMAP-Sample-Clusters.pdf", ArchRProj = subProj, addDOC = FALSE, width = 5, height = 5)

# Save unfiltered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Remove doublet clusters
##########################################################################################

# There are a few clusters that appear to be mostly doublets / poor quality cells
# (Poor TSS enrichment, few marker peaks / GS in downstream analysis, higher cellCallUncertainty)
# and were not filtered by automated doublet removal or by basic QC filtering
# They will be manually filtered here.
proj <- loadArchRProject(paste0(wd, "/naive_filtered_output/"), force=TRUE)

# Identify Marker Gene through Pairwise Test vs Bias-Matched Background:
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

# Marker genes we want to highlight for labeling broad clusters:
markerGenes  <- c(
  "CD34", "HLF", "AVP", "SPINK2","AZU1","ELANE",  # HSPC
  "CLC","GZMB",
  "CEBPE","ITGAM", "MMP9", "PGLYRP1", "ORM1", "CEACAM8","PADI4",  # Neutrophils
  "CD3D", "CD8A", "CD4", "PTPRC", "FOXP3", "IKZF2", "CCL5","GZMK", # T-cells
  "NKG7","GNLY", #NK
  "HBA2","GYPA","AHSP",#Ery
  "CD19", "MS4A1","PAX5","JCHAIN","RAG1", # B-cells
  "ITGAX", "CD1C", "CD1A", "CLEC1A", "CD207", # Dendritic cells (ITGAX = cd11c)
  "CD14", "CD86", "CD74", "CD163" #Monocytes / macrophages
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(hm, name = "filtered-GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

#visualize doublet score across all clusters
plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "cellColData", 
  name = "DoubletScore",
  plotAs = "violin",
  alpha = 0.4
)

plotEmbedding(proj, colorBy="cellColData", name="DoubletScore", 
              embedding = "mnnUMAP", plotAs="points", size=pointSize, 
              labelMeans=FALSE, imputeWeights=getImputeWeights(proj))

# Remove clusters that have poor quality (enriched for high doublet score cells, poor cell quality, incompatible marker gene scores, etc.)
nonMultipletCells <- getCellNames(proj)[proj$Clusters %ni% c("C7")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = nonMultipletCells,
  outputDirectory = "multiplets_removed_output",
  dropCells=TRUE, force=TRUE
)
saveArchRProject(proj)

# Now, redo clustering and visualization:

# Reduce Dimensions with Iterative LSI
set.seed(1)
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI",
    sampleCellsPre = 20000,
    varFeatures = 50000, 
    dimsToUse = 1:50,
    force = TRUE
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  dimsToUse = 2:30,
  name = "Harmony",
  groupBy = "Sample",
  force = T
)

fMNN.res <- mnnCorrect(t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'HC_1'),1:50]),
                       t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'HC_2'),1:50]),
                       t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'ST_4'),1:50]),
                       t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'ST_5'),1:50]),
                       t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'ST_6'),1:50]),
                       t(proj@reducedDims@listData[["IterativeLSI"]]@listData[["matSVD"]][which(proj$Sample == 'ST_8'),1:50]))

proj@reducedDims[['LSI_MNN']] <- SimpleList(
  matDR = t(fMNN.res@assays@data$corrected), 
  params = NA,
  date = Sys.time(),
  scaleDims = NA, #Do not scale dims after
  corToDepth = NA
)

# Identify Clusters from Iterative LSI/LSI_MNN
proj <- addClusters(
    input = proj,
    reducedDims = "LSI_MNN",
    dimsToUse = 2:30,
    method = "Seurat",
    name = "Clusters",
    maxClusters = 30,
    resolution = 0.8,
    force = TRUE
)

set.seed(1)
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 60, 
    dimsToUse = 2:30,
    minDist = 0.6, 
    metric = "cosine",
    force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "HarmonyUMAP", 
  nNeighbors = 60, 
  dimsToUse = 2:30,
  minDist = 0.6, 
  metric = "cosine",
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "LSI_MNN", 
  name = "mnnUMAP", 
  nNeighbors = 60, 
  dimsToUse = 2:30,
  minDist = 0.6, 
  metric = "cosine",
  force = TRUE
)

samp_cmap <- c("#FF7D7D","#ba79b1","#896d47","#86908a","#e18a3a","#a45e44")
names(samp_cmap) <- names(table(proj$Sample))
disease_cmap <- c("#d3a237","#4f794a")
names(disease_cmap) <- names(table(proj$group_info))
# Relabel clusters so they are sorted by cluster size
proj <- relabelClusters(proj)
proj <- addImputeWeights(proj)
# Make various cluster plots:
proj <- visualizeClustering(proj, pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)
#don't forget to remove repeated barplots which are the same as above
proj <- visualizeClustering(proj, prefix = 'mnn', embedding = 'mnnUMAP',
                            pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "HarmonyUMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "HarmonyUMAP")

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "mnnUMAP")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "mnnUMAP")

plotPDF(p1,p2,p3,p4, name = "Plot-BatchCorrection-UMAP-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Save filtered ArchR project
saveArchRProject(proj)

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# All valid gene names
geneNames <- rowData(markersGS)$name

# Lists of 'marker genes'
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.00")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.00", 
  labelMarkers = markerGenes,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(hm, name = "GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

# Marker gene imputation with Magic:
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "mnnUMAP",
  imputeWeights = getImputeWeights(proj), 
  plotAs="points", size = pointSize
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  tileSize=250,
  minCells=25
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

bmrna <- readRDS('~/Desktop/ATAC/bmmc_sletp_rna.rds')

####not yet
library(zellkonverter)
all_h5ad <- readH5AD('../../bmmc_annot.h5ad',reader = "python")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "LSI_MNN",
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  seRNA = all_h5ad,
  addToArrow = FALSE,
  groupRNA = "cell_type_eryMK",
  nameCell = "predictedCell_rna",
  nameGroup = "predictedGroup_rna",
  nameScore = "predictedScore_rna"
)
pal <- paletteDiscrete(values = all_h5ad$cell_type_eryMK)

p_pred <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup_rna", 
  embedding = "mnnUMAP",
  pal = pal
)

#####

plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = "HLA-DPA1", 
  embedding = "mnnUMAP",
  imputeWeights = getImputeWeights(proj), 
  plotAs="points", size = pointSize
)
# Assign cluster names based on marker genes:
clustNames <- list(
  "C1" = "T",
  "C2" = "Myeloid", 
  "C3" = "T", 
  "C4" = "Neutrophil",
  "C5" = "T", 
  "C6" = "Myeloid",
  "C7" = "Neutrophil",
  "C8" = "T",
  "C9" = "T",
  "C10" = "NK",
  "C11" = "naive_B",
  "C12" = "Myeloid",
  "C13" = "Erythroid",
  "C14" = "memory_B",
  "C15" = "MEP",
  "C16" = "T",
  "C17" = "Plasma",
  "C18" = "Myeloid",
  "C19" = "HSC_MPP",
  "C20" = "Myeloid",
  "C21" = "NeuP",
  "C22" = "preB",
  "C23" = "EBMP",
  "C24" = "proB",
  "C25" = "NK",
  "C26" = "Myeloid"
)
proj$NamedClust <- clustNames[proj$Clusters] %>% unlist()
p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "NamedClust", embedding = "mnnUMAP")
plotPDF(p, name = "Plot-majorcelltype-UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

plotDir <- paste0(proj@projectMetadata$outputDirectory, "/Plots")
barwidth=0.9
namedClustBySamp <- fractionXbyY(proj$NamedClust, proj$Sample, add_total=TRUE, xname="NamedClust", yname="Sample")
pdf(paste0(plotDir, "/sampleByNamedClustBarPlot.pdf"))
print(stackedBarPlot(namedClustBySamp, cmap=samp_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

namedClustByGroup <- fractionXbyY(proj$NamedClust, proj$group_info, add_total=TRUE, xname="NamedClust", yname="Sample")
pdf(paste0(plotDir, "/groupByNamedClustBarPlot.pdf"))
print(stackedBarPlot(namedClustByGroup, cmap=disease_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()
saveArchRProject(proj)

##########################################################################################
# Subcluster all major cell groups
##########################################################################################

# All major cell clusters will be subclustered for downstream analysis
# Peaks will also be called on subclustered groups

subClusterGroups <- list(
  "NK_T" = c("NK","T"), 
  "HSPC" = c("Erythroid","HSC_MPP","EBMP","MEP","NeuP"),
  "Myeloid" = c("Myeloid"),
  "Neutrophil" = c("Neutrophil"),
  "B_lineage" = c("proB","preB","naive_B","memory_B","Plasma")
)

subClusterCells <- lapply(subClusterGroups, function(x){
  getCellNames(proj)[as.character(proj@cellColData[["NamedClust"]]) %in% x]
})

subClusterArchR <- function(proj, subCells, outdir){
  # Subset an ArchR project for focused analysis
  
  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
    ArchRProj = proj,
    cells = subCells,
    outputDirectory = outdir,
    dropCells = TRUE,
    force = TRUE
  )
  saveArchRProject(sub_proj)
  
  # Delete stuff we don't want to copy...
  unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
  unlink(paste0(outdir, "/IterativeLSI/*"), recursive = TRUE)
  unlink(paste0(outdir, "/Embeddings/*"), recursive = TRUE)
  unlink(paste0(outdir, "/ImputeWeights/*"), recursive = TRUE)
}

subgroups <- names(subClusterGroups)

sub_proj_list <- lapply(subgroups, function(sg){
  message(sprintf("Subsetting %s...", sg))
  outdir <- sprintf("~/Desktop/ATAC/202312/subclustered_%s", sg)
  subClusterArchR(proj, subCells=subClusterCells[[sg]], outdir=outdir)
})
names(sub_proj_list) <- subgroups

#####
#more detailed cell sub-type will be added to proj after re-clustering in each lineage

motif_matrix <- getMatrixFromProject(atac_proj,useMatrix = 'MotifMatrix')
matrix <- motif_matrix@assays@data$z
TF <- c('PAX5','EBF1','ID3','ETS2','FLI1','GATA1','GATA2',"PBX1","SPI1",'LMO4','IRF7','IRF8')
TF_num <- c()
for (i in TF){
  TF_num <- c(TF_num,grep(i,rownames(matrix)))
}
TF_num <- setdiff(TF_num,grep('SREBF1',rownames(matrix)))

atac_proj$group_celltype <- paste0(atac_proj$NamedClust,'(',atac_proj$group_info,')')

tf_matrix <- matrix(nrow = length(TF_num), ncol = length(table(atac_proj$group_celltype)))
rownames(tf_matrix) <- motif_name
colnames(tf_matrix) <- names(table(atac_proj$group_celltype))
motif_name <- rownames(matrix)[TF_num]

for (i in names(table(atac_proj$group_celltype))){
  for (j in motif_name){
    cells <- getCellNames(atac_proj)[which(atac_proj$group_celltype == i)]
    tf_matrix[j,i] <- mean(matrix[j,cells])
  }
}

genescore_matrix <- getMatrixFromProject(atac_proj,useMatrix = 'GeneScoreMatrix')
matrix2 <- genescore_matrix@assays@data$GeneScoreMatrix
gene_names <- genescore_matrix@elementMetadata@listData[["name"]]
rownames(matrix2) <- gene_names
TF_num2 <- c()
for (i in TF){
  TF_num2 <- c(TF_num2,grep(i,rownames(matrix2)))
}
gene_name <- intersect(rownames(matrix2)[TF_num2],TF)

tf_matrix2 <- matrix(nrow = length(gene_name), ncol = length(table(atac_proj$group_celltype)))
rownames(tf_matrix2) <- gene_name
colnames(tf_matrix2) <- names(table(atac_proj$group_celltype))

for (i in names(table(atac_proj$group_celltype))){
  for (j in gene_name){
    cells <- getCellNames(atac_proj)[which(atac_proj$group_celltype == i)]
    tf_matrix2[j,i] <- mean(matrix2[j,cells])
  }
}
library(pheatmap)
p1 <- pheatmap::pheatmap(tf_matrix,scale = 'row',cellwidth = 10,cellheight = 10,
                         cluster_rows = F,cluster_cols = F,
                         main = 'Transcription Factor Motif Activity')

p2 <- pheatmap::pheatmap(tf_matrix2,scale = 'row',cellwidth = 10,cellheight = 10,
                         cluster_rows = F,cluster_cols = F,
                         main = 'Transcription Factor Gene Score')

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")
pdf(paste0(plotDir, "/TF-certain-Motif-heatmap.pdf"),width = 10,height = 10)
p1
dev.off()
pdf(paste0(plotDir, "/TF-certain-GeneScore-heatmap.pdf"),width = 10,height = 10)
p2
dev.off()


