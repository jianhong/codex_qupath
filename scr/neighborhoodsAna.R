#!/usr/bin/env Rscript
# This script allows to detect the cellstype from the output of detectCellType.R
#
# Script author: Jianhong Ou
# Script should work after running detectCellType.R
# This source code is licensed under the MIT license
# How to run in command line:
# Rscript neighborhoodsAna path/to/seuratObject.rds
# How to run in Rstudio/R:

args = commandArgs(trailingOnly=TRUE)
# args = c(path='path/to/seuratObject.rds')
# eg args <- '/Users/ouj/Downloads/qupath/run_20230926/measurementsExport13635325425080860669/seuratObj/Julia Visguass_RealRun(TMA1)_230731.preprocessed.cell.median.seu.withCellType.rds'

size <- 1 # plot dot size
alpha <- 0.8 # alpha of the dot

library(Seurat)
library(ggplot2)
library(hoodscanR)
library(SpatialExperiment)
library(scico)
library(ComplexHeatmap)

pilot <- readRDS(args[1])
pf <- dirname(args[1])
spe <- SpatialExperiment(
    assay = GetAssayData(pilot, assay = 'GMM', layer = 'data'),
    colData = pilot[[]],
    spatialCoordsNames = c('Cell.X', 'Cell.Y')
)
p <- plotTissue(spe, color = celltype, size = size, alpha = alpha)
ggsave(file.path(pf, 'celltypePlot.pdf'), plot = p, 
       width = 14, height = 14)

# neighborhoods distance scanning, k = maximal neighbor cells
fnc <- findNearCells(spe, k = 25, anno_col = 'celltype')
# probability matrix
pm <- scanHoods(fnc$distance)
hoods <- mergeByGroup(pm, fnc$cells)
spe <- mergeHoodSpe(spe, hoods)
spe <- calcMetrics(spe, pm_cols = colnames(hoods))
p <- plotTissue(spe, color = entropy, size = size, alpha = alpha) +
    scale_color_scico(palette = 'berlin')
ggsave(file.path(pf, 'entropy.of.mixture.neighborhood.pdf'), plot = p, 
       width = 14, height = 14)
p <- plotTissue(spe, color = perplexity, size = size, alpha = alpha) +
    scale_color_scico(palette = 'berlin')
ggsave(file.path(pf, 'perplexity.of.mixture.neighborhood.pdf'), plot = p, 
       width = 14, height = 14)
p <- plotColocal(spe, pm_cols = colnames(hoods))
pdf(file.path(pf, 'neighborhoodDistancePearsonCoor.pdf'))
draw(p) 
dev.off()
clusterK <- 12
spe <- clustByHood(spe, pm_cols = colnames(hoods), k=clusterK)
p <- plotProbDist(spe, pm_cols = colnames(hoods), 
                  by_cluster = TRUE, plot_all = TRUE, 
                  show_clusters = as.character(seq(clusterK)))
ggsave(file.path(pf, 'ProbDistClust.neighborhood.pdf'), plot = p, 
       width = 14, height = 14)
p <- plotTissue(spe, color = clusters, size = size, alpha = alpha)
ggsave(file.path(pf, 'celltypeProbDistClust.pdf'), plot = p, 
       width = 14, height = 14)
saveRDS(spe, file.path(pf, 'neighborhood.spe.rds'))
