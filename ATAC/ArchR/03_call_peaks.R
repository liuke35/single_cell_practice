suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
})

source('~/Desktop/ATAC/analysis_function/archr_helper.R')
source('~/Desktop/ATAC/analysis_function/plotting_config.R')
source('~/Desktop/ATAC/analysis_function/misc_helper.R')
source('~/Desktop/ATAC/analysis_function/matrix_helper.R')
source('~/Desktop/ATAC/analysis_function/homemade_archr.R')

addArchRThreads(threads = 4)

wd <- "~/Desktop/ATAC/202312/01_data_integration/multiplets_removed_output/"
outdir <- "~/Desktop/ATAC/202312/fine_clustered"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
pointSize <- 0.25

orig_proj <- loadArchRProject(wd, force=TRUE)

atac_proj <- subsetArchRProject(
  ArchRProj = orig_proj,
  cells = getCellNames(orig_proj),
  outputDirectory = outdir,
  dropCells = TRUE,
  force = TRUE
)
saveArchRProject(atac_proj)
rm(orig_proj); gc()

# Delete stuff we don't need to duplicate
unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
unlink(paste0(outdir, "/GroupCoverages"), recursive = TRUE)
unlink(paste0(outdir, "/PeakCalls"), recursive = TRUE)
unlink(paste0(outdir, "/*_filtered_barcodes.txt"), recursive = TRUE)

##########################################################################################
# Identify cluster labels from subclustered ArchR projects
##########################################################################################

atac_proj <- loadArchRProject(outdir, force=TRUE)
# color palettes
sample_cmap = readRDS('~/Desktop/ATAC/202312/01_data_integration/samp_cmap.rds')
disease_cmap <- readRDS('~/Desktop/ATAC/202312/01_data_integration/disease_cmap.rds')

subclustered_projects <- c("NK_T", "Myeloid", "HSPC", "Neutrophil", "B_lineage")

FineClustLabels <- atac_proj$NamedClust # Default to NamedClust where no subgroup label exists
names(FineClustLabels) <- getCellNames(atac_proj)

for(subgroup in subclustered_projects){
  message(sprintf("Reading in subcluster %s", subgroup))
  # Read in subclustered project
  sub_dir <- sprintf("~/Desktop/ATAC/202312/subclustered_%s", subgroup)
  sub_proj <- loadArchRProject(sub_dir, force=TRUE)
  
  # Add FineClust to full ArchR project
  FineClustLabels[getCellNames(sub_proj)] <- sub_proj$FineClust
}

# Add FineCluster information to full ArchR project
atac_proj$FineClust <- FineClustLabels[getCellNames(atac_proj)]

p1 <- plotEmbedding(ArchRProj=atac_proj, colorBy = "cellColData", name="NamedClust",
                    embedding = "mnnUMAP", plotAs="points", size=pointSize, labelMeans = FALSE)
p2 <- plotEmbedding(ArchRProj = atac_proj, colorBy="cellColData", name="FineClust",
                    embedding="mnnUMAP", plotAs="points", size=pointSize, labelMeans=FALSE)
ggAlignPlots(p1,p2, type="h")
plotPDF(p1,p2, name = "Plot-UMAP-FineClusters.pdf", ArchRProj=atac_proj, addDOC=FALSE, width=5, height=5)

saveArchRProject(atac_proj)

##########################################################################################
# Call Peaks
##########################################################################################
library(BSgenome.Hsapiens.UCSC.hg38)

atac_proj <- loadArchRProject(outdir, force=TRUE)

# Create Group Coverage Files that can be used for downstream analysis
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineClust", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)

# Find Path to Macs2 binary
pathToMacs2 <- '/Users/keliu/anaconda3/bin/macs2'

# Call Reproducible Peaks w/ Macs2
atac_proj <- addReproduciblePeakSet(
  ArchRProj = atac_proj, 
  groupBy = "FineClust", 
  peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

# Add Peak Matrix
atac_proj <- addPeakMatrix(ArchRProj = atac_proj, force = TRUE)

# Save project
saveArchRProject(atac_proj)

