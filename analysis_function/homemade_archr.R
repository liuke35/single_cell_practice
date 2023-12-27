suppressPackageStartupMessages({
  library(ArchR)
  library(cowplot)
  library(ggplot2)
})

pseudotime_plot <- function(proj, trajectory, group, features, useMatrix = 'GeneScoreMatrix'){
  cell <- which(is.na(getCellColData(proj,trajectory))== F)
  matrix <- getMatrixFromProject(proj,useMatrix = useMatrix)
  if (useMatrix == 'GeneScoreMatrix'){
    gene_names <- matrix@elementMetadata@listData[["name"]]
    matrix <- matrix@assays@data$GeneScoreMatrix
    rownames(matrix) <- gene_names
  }else if(useMatrix == 'MotifMatrix'){
    matrix <- matrix@assays@data$z
  }
  df <- data.frame(cellname = getCellNames(proj)[cell],
                   pseudotime = getCellColData(proj,trajectory)[,1][cell],
                   group_info = getCellColData(proj,group)[,1][cell])
  for (i in features){
    df[i] <- matrix[i,df$cellname]
  }
  
  gg_list <- list()
  for (i in 1:length(features)){
    feat <- data.frame(pseudotime = df$pseudotime,group_info = df$group_info,
                       gene = df[plot_features[i]])
    colnames(feat) <- c('pseudotime','group_info','gene')
    p <- ggplot(feat, aes(x = pseudotime, y = gene, color = group_info)) +
      geom_point(size = .5) +
      geom_smooth(method = 'loess', se = F)+
      xlab("Pseudotime") + 
      ylab(features[i])+
      theme(panel.grid = element_blank(),
            plot.background = element_rect(fill = 'white'),
            panel.background = element_rect(fill = 'white'),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            legend.position = 'right')
    gg_list[[i]] <- p
  }
  return(gg_list)
}

motif_plot <- function(atac_proj,markerTest,prefix){
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  
  ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = atac_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
  
  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  
  ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  plotPDF(ggUp, ggDo, name = paste0(prefix,"-Markers-Motifs-Enriched"), width = 5, height = 5, ArchRProj = atac_proj, addDOC = FALSE)
}
