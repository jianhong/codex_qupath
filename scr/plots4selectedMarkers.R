#!/usr/bin/env Rscript
# This script allows to detect the cells and merge multinucleis annotations
#
# Script author: Jianhong Ou
# Script should work after running createSeuratObj.R
# This source code is licensed under the MIT license
# How to run in command line:
# Rscript createSeuratObj.R path/to/preprocessed/seu/rds [selected markers ...]

args = commandArgs(trailingOnly=TRUE)
# How to run in Rstudio/R:
# uncomment the following line and change the preprocessed/seu/rds
# args <- c('~/Downloads/qupath/run_20230926/measurementsExport13635325425080860669/seuratObj/Julia Visguass_RealRun(TMA1)_230731.preprocessed.cell.median.seu.rds')

library(Seurat)
library(ggplot2)

pilot <- readRDS(args[1])

if(length(args)>1){
    selected.genes <- args[-1]
}else{
    selected.genes <- VariableFeatures(pilot)
}
## if you want to change the selected markers, uncomment the following line and change the markers
# selected.genes <- c('"cd34", "lcp1", "hif1a", "cd45")

pt.size <- ifelse(ncol(pilot)>500, 0, .5)
pf <- dirname(args[1])
prefix <- sub(".seu.rds", '', basename(args[1]))
p <- DotPlot(pilot, features = rev(selected.genes), idents=pilot$seurat_clusters) + coord_flip()
ggsave(file.path(pf, paste0(prefix, '.selected.genes.dotplot.pdf')), plot=p, width=12)

for(layer in Layers(pilot)){
    cor <- cor(t(as.matrix(GetAssayData(pilot, assay = 'RNA', layer = layer))[selected.genes, ]), method='spearman')
    pdf(file.path(pf, paste0(prefix, '.selected.genes.correlation.by.', layer, '.pdf')))
    p <- heatmap(cor, scale = 'none')
    dev.off()
}
