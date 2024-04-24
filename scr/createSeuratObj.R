#!/usr/bin/env Rscript
# This script allows to detect the cells and merge multinucleis annotations
#
# Script author: Jianhong Ou
# Script should work after running ExportCellDetectionMeasurement.groovy
# This source code is licensed under the MIT license
# How to run in command line:
# Rscript createSeuratObj.R path/to/measurementsExport/folder
# How to run in Rstudio/R:
# uncomment line 18 and change the folder to the measurements folder
# for normalization, available method are
# 'LogNormalize', 'CLR' (centered log ratio), 'RC' (CPM), 'zscore'  
# Assign nFeatures as a number will select top N features with maximal variance
# Assign avoid.features as a character vector by the markers will decide to avoid those markers as features for PCA analysis
# eg avoid.features = c('sox9', 'hist3h2a')

args = commandArgs(trailingOnly=TRUE)
# args = c(path='path/to/measurementsExport/folder', useValue='Median', normalizationMethod='zscore', nFeatures = 10, avoid.features=NULL)

library(Seurat)
library(ggplot2)
library(dittoSeq)

fs <- dir(args[1], "qptiff.tsv", full.names=TRUE)
useValue <- ifelse(length(args)>1, args[2], 'Median')
normalizationMethod <- ifelse(length(args)>2, args[3], 'zscore')
nFeatures <- ifelse(length(args)>3, as.numeric(args[4]), 10)
#asyN <- c('Nucleus', 'Membrane', 'Cytoplasm', 'Cell')
#asyN2 <- c('Mean', 'Median', 'Min', 'Max', 'Std.Dev')
markerLocations <- c(
    'DAPI.01'='Nucleus',
    'CD163'='Cell',
    'CD127'='Cell',
    'sox9'='Nucleus',
    'CD14'='Cell',
    'CD4'='Cell',
    'mpo'='Cell',
    'cd68'='Cell',
    'cd45'='Cell',
    'foxp3'='Cell',
    'lcp1'='Cell',
    'pd1'='Cell',
    'cd8'='Cell',
    'hif1a'='Cell',
    'ki67'='Nucleus',
    'cd19'='Cell',
    'cd206'='Cell',
    'znf90'='Cell',
    'asma'='Cell',
    'cd3'='Cell',
    'cd86'='Cell',
    'hist3h2a'='Cell',
    'cd34'='Cell',
    'CD11c'='Cell'
)
avoid.features <- NULL
if(length(args)>4){
    avoid.features <- args[5]
    if(length(avoid.features)){
        stopifnot(all(avoid.features %in% names(markerLocations)))
    }
}


## meta data information filter
metaInfoStartWith <- '^(Nucleus|Cell)'

for(f in fs){
    message('Working on:', f)
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
    dx$seurat_obj_raw <- dx[[paste0("Cell_", useValue)]]
    for(marker in names(markerLocations)){
        dx$seurat_obj_raw[marker, ] <- 
            dx[[paste(markerLocations[marker], useValue, sep="_")]][marker, ]
    }
    pff <- file.path(pf, 'counts')
    dir.create(pff)
    null <- mapply(dx, names(dx), FUN=function(.dat, .name){
        write.csv(.dat, file.path(pff, paste0(prefix, '.', .name, '.csv')))
    })
    seu <- CreateSeuratObject(dx$seurat_obj_raw,
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
    
    if(normalizationMethod=='zscore'){
        new.data <- scale(GetAssayData(seu, assay = 'RNA', layer = 'counts'))
        pilot <- SetAssayData(seu, assay = 'RNA', layer='data', new.data = new.data)
    }else{
        pilot <- NormalizeData(seu, normalization.method = normalizationMethod)
    }
    pilot <- FindVariableFeatures(pilot, selection.method = "vst", nfeatures = nFeatures)
    if(length(avoid.features)){
        pilot$RNA@meta.data$var.features.rank[pilot$RNA@meta.data$var.features %in% avoid.features] <- NA
        pilot$RNA@meta.data$var.features[pilot$RNA@meta.data$var.features %in% avoid.features] <- NA
    }
    plot1 <- VariableFeaturePlot(pilot)
    plot2 <- LabelPoints(plot = plot1, 
                         points = VariableFeatures(pilot),
                         repel = TRUE)
    ggsave(file.path(pf, paste0(prefix, 'variableFeatures.pdf')), plot=plot2)
    
    all.genes <- rownames(pilot)
    pilot <- ScaleData(pilot, features = all.genes)
    pilot <- RunPCA(pilot, features = VariableFeatures(pilot))
    dim <- seq.int(ncol(Embeddings(pilot$pca)))
    pilot <- FindNeighbors(pilot, dims = dim)
    pilot <- FindClusters(pilot, resolution = 0.3)
    pilot <- RunUMAP(pilot, reduction = "pca", dims = dim)
    saveRDS(pilot, file.path(pf, paste0(prefix, ".preprocessed.cell.", useValue,".seu.rds")))
    pt.size <- ifelse(ncol(pilot)>500, 0, .5)
    p <- VlnPlot(pilot, features = c("nFeature_RNA", "nCount_RNA", all.genes), pt.size = pt.size)
    ggsave(file.path(pf, paste0(prefix, '.vlnplot.pdf')), plot=p, width=28, height = 28)
    p <- DotPlot(pilot, features = rev(all.genes), idents=pilot$seurat_clusters) + coord_flip()
    ggsave(file.path(pf, paste0(prefix, '.dotplot.pdf')), plot=p, width=12)
    p <- DimPlot(pilot, label = TRUE)
    ggsave(file.path(pf, paste0(prefix, '.dimplot.pdf')), plot=p, width=9)
    p <- DimPlot(pilot, reduction = 'pos', label = TRUE) + scale_y_reverse()
    ggsave(file.path(pf, paste0(prefix, '.positionPlot.pdf')), plot=p, width=9)
    ncol <- ceiling(sqrt(length(unique(pilot$seurat_clusters))))
    p <- DimPlot(pilot, reduction = 'pos', label = FALSE, split.by = 'seurat_clusters', ncol = ncol, pt.size = 2) + scale_y_reverse()
    ggsave(file.path(pf, paste0(prefix, '.positionPlot.by.cluster.pdf')), plot=p, width=6*ncol, height = 6*ncol)
    p <- dittoBarPlot(pilot, var = 'Cell.isMultiNucleis', group.by='seurat_clusters')
    ggsave(file.path(pf, paste0(prefix, '.multiNucleis.vs.clusters.pdf')), plot=p, width=7)
    p <- DimPlot(pilot, reduction = 'pos', label = TRUE, split.by = 'Cell.isMultiNucleis', pt.size = 2) + scale_y_reverse()
    ggsave(file.path(pf, paste0(prefix, '.positionPlot.by.isMultiNucleis.pdf')), plot=p, width=14, height = 7)
    p <- VlnPlot(pilot, features = colnames(meta)[vapply(meta, is.numeric, logical(1L))], pt.size=pt.size)
    ggsave(file.path(pf, paste0(prefix, '.metadata.by.cluster.pdf')), plot=p, width=28, height = 28)
    p <- VlnPlot(pilot, features = colnames(meta)[vapply(meta, is.numeric, logical(1L))], group.by = 'Cell.isMultiNucleis', pt.size = pt.size)
    ggsave(file.path(pf, paste0(prefix, '.metadata.by.classification.pdf')), plot=p, width=28, height = 28)
    for(layer in Layers(pilot)){
        cor <- cor(t(as.matrix(GetAssayData(pilot, assay = 'RNA', layer = layer))), method='spearman')
        pdf(file.path(pf, paste0(prefix, '.correlation.by.', layer, '.pdf')))
        p <- heatmap(cor, scale = 'none')
        dev.off()
    }
    expMeta <- data.frame('Object ID'=colnames(pilot), 'Class'=paste0('Cluster', pilot$seurat_clusters))
    write.csv(expMeta, file.path(pf, paste0(prefix, '.seurat_cluster.scv')), quote = FALSE, row.names = FALSE)
}

