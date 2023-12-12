#!/usr/bin/env Rscript
# This script allows to detect the cells and merge multinucleis annotations
#
# Script author: Jianhong Ou
# Script should work after running ExportCellDetectionMeasurement.groovy
# This source code is licensed under the MIT license
# Howto run in command line:
# Rscript createSeuratObj.R path/to/measurementsExport/folder

args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(ggplot2)
library(dittoSeq)

fs <- dir(args[1], "qptiff.tsv", full.names=TRUE)
#asyN <- c('Nucleus', 'Membrane', 'Cytoplasm', 'Cell')
#asyN2 <- c('Mean', 'Median', 'Min', 'Max', 'Std.Dev')
metaInfoStartWith <- '^(Nucleus|Cell)'

for(f in fs){
    pf <- file.path(dirname(f), 'seuratObj')
    prefix <- sub(".qptiff.tsv", '', basename(f))
    dir.create(pf)
    d <- read.delim(f)
    meta <- d[, grepl(metaInfoStartWith, colnames(d))]
    rownames(meta) <- meta$Cell.ID
    meta$Cell.ID <- NULL
    meta$Cell.isMultiNucleis <- ifelse(grepl('MultiNucleis', meta$Cell.classification), 'MultiNucleis', 'SingleNuclei')
    d <- d[, !grepl(metaInfoStartWith, colnames(d))]
    cn <- colnames(d)
    genes <- sub("..Cell..Median", '', cn[grepl('Cell..Median', cn)])
    cn1 <- unique(sub('^(.*?)\\.\\.(Nucleus|Membrane|Cytoplasm|Cell)\\.\\.(.*?)$', "..\\2..\\3", cn))
    dx <- lapply(cn1, function(.ele){
        dat <- d[, cn[grepl(.ele, cn)]]
        colnames(dat) <- sub(.ele, '', colnames(dat))
        dat <- t(dat)
        colnames(dat) <- rownames(meta)
        dat
    })
    names(dx) <- sub("^_", "", gsub("\\.\\.", '_', cn1))
    pff <- file.path(pf, 'counts')
    dir.create(pff)
    null <- mapply(dx, names(dx), FUN=function(.dat, .name){
        write.csv(.dat, file.path(pff, paste0(prefix, '.', .name, '.csv')))
    })
    seu <- CreateSeuratObject(dx[["Cell_Median"]],
                              project = prefix,
                              assay = "RNA",
                              meta.data = meta,
                              min.cells = 3,
                              min.features = 1)
    pos <- meta[, c('Cell.X', 'Cell.Y')]
    colnames(pos) <- paste('POS_', 1:2)
    seu@reductions$pos <- CreateDimReducObject(embeddings = as(pos, "matrix"),
                                key = 'POS_',
                                assay = DefaultAssay(seu))
    
    saveRDS(seu, file.path(pf, paste0(prefix, ".origin.cell.median.seu.rds")))

    pilot <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    pilot <- FindVariableFeatures(pilot, selection.method = "vst")
    all.genes <- rownames(pilot)
    pilot <- ScaleData(pilot, features = all.genes)
    pilot <- RunPCA(pilot, npcs = 15)
    pilot <- FindNeighbors(pilot, dims = 1:15)
    pilot <- FindClusters(pilot, resolution = 0.3)
    pilot <- RunUMAP(pilot, reduction = "pca", dims = 1:15)
    saveRDS(pilot, file.path(pf, paste0(prefix, ".preprocessed.cell.median.seu.rds")))
    p <- VlnPlot(pilot, features = c("nFeature_RNA", "nCount_RNA", all.genes))
    ggsave(file.path(pf, paste0(prefix, '.vlnplot.pdf')), plot=p, width=28, height = 28)
    p <- DotPlot(pilot, features = rev(all.genes), idents=pilot$seurat_clusters) + coord_flip()
    ggsave(file.path(pf, paste0(prefix, '.dotplot.pdf')), plot=p, width=12)
    p <- DimPlot(pilot, label = TRUE)
    ggsave(file.path(pf, paste0(prefix, '.dimplot.pdf')), plot=p, width=9)
    p <- DimPlot(pilot, reduction = 'pos', label = TRUE)
    ggsave(file.path(pf, paste0(prefix, '.positionPlot.pdf')), plot=p, width=9)
    p <- dittoBarPlot(pilot, var = 'Cell.isMultiNucleis', group.by='seurat_clusters')
    ggsave(file.path(pf, paste0(prefix, '.multiNucleis.vs.clusters.pdf')), plot=p, width=7)
    p <- VlnPlot(pilot, features = colnames(meta)[vapply(meta, is.numeric, logical(1L))])
    ggsave(file.path(pf, paste0(prefix, '.metadata.by.cluster.pdf')), plot=p, width=28, height = 28)
    p <- VlnPlot(pilot, features = colnames(meta)[vapply(meta, is.numeric, logical(1L))], group.by = 'Cell.isMultiNucleis')
    ggsave(file.path(pf, paste0(prefix, '.metadata.by.classification.pdf')), plot=p, width=28, height = 28)
}

