#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(ggplot2)
fs <- dir(args[1], "qptiff.tsv", full.names=TRUE)
#asyN <- c('Nucleus', 'Membrane', 'Cytoplasm', 'Cell')
#asyN2 <- c('Mean', 'Median', 'Min', 'Max', 'Std.Dev')

for(f in fs){
    pf <- file.path(dirname(f), 'seuratObj')
    prefix <- sub(".qptiff.tsv", '', basename(f))
    dir.create(pf)
    d <- read.delim(f)
    meta <- d[, 1:13]
    cn <- colnames(d)
    genes <- sub("..Cell..Median", '', cn[grepl('Cell..Median', cn)])
    cn1 <- cn[-c(1:13)]
    cn1 <- unique(sub('^(.*?)\\.\\.(Nucleus|Membrane|Cytoplasm|Cell)\\.\\.(.*?)$', "..\\2..\\3", cn1))
    dx <- lapply(cn1, function(.ele){
        dat <- d[, cn[grepl(.ele, cn)]]
        colnames(dat) <- sub(.ele, '', colnames(dat))
        dat <- t(dat)
        colnames(dat) <- paste0('cell_', seq.int(ncol(dat)))
        dat
    })
    rownames(meta) <- paste0('cell_', seq.int(nrow(meta)))
    names(dx) <- sub("^_", "", gsub("\\.\\.", '_', cn1))
    pff <- file.path(pf, 'counts')
    dir.create(pff)
    mapply(dx, names(dx), FUN=function(.dat, .name){
        write.csv(.dat, file.path(pff, paste0(prefix, '.', .name, '.csv')))
    })
    seu <- CreateSeuratObject(dx[["Cell_Median"]],
                              project = prefix,
                              assay = "RNA",
                              meta.data = meta,
                              min.cells = 3,
                              min.features = 3)
    saveRDS(seu, file.path(pf, paste0(prefix, ".origin.seu.rds")))
    
    pilot <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    pilot <- FindVariableFeatures(pilot, selection.method = "vst")
    all.genes <- rownames(pilot)
    pilot <- ScaleData(pilot, features = all.genes)
    pilot <- RunPCA(pilot, npcs = 15)
    pilot <- FindNeighbors(pilot, dims = 1:15)
    pilot <- FindClusters(pilot, resolution = 0.3)
    pilot <- RunUMAP(pilot, reduction = "pca", dims = 1:15)
    saveRDS(pilot, file.path(pf, paste0(prefix, ".preprocessed.seu.rds")))
    p <- VlnPlot(pilot, features = c("nFeature_RNA", "nCount_RNA", all.genes))
    ggsave(file.path(pf, paste0(prefix, '.vlnplot.pdf')), plot=p, width=28)
    p <- DotPlot(pilot, features = rev(all.genes), idents=pilot$seurat_clusters) + coord_flip()
    ggsave(file.path(pf, paste0(prefix, '.dotplot.pdf')), plot=p, width=12)
    p <- DimPlot(pilot, label = TRUE)
    ggsave(file.path(pf, paste0(prefix, '.dimplot.pdf')), plot=p, width=9)
}