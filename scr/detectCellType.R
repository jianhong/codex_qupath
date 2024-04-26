#!/usr/bin/env Rscript
# This script allows to detect the cellstype from the output of createSeuratObj.R
#
# Script author: Jianhong Ou
# Script should work after running createSeuratObj.R
# This source code is licensed under the MIT license
# How to run in command line:
# Rscript detectCellType.R path/to/seuratObject.rds
# How to run in Rstudio/R:

args = commandArgs(trailingOnly=TRUE)
# args = c(path='path/to/seuratObject.rds')
# eg args <- '/Users/ouj/Downloads/qupath/run_20230926/measurementsExport13635325425080860669/seuratObj/Julia Visguass_RealRun(TMA1)_230731.preprocessed.cell.median.seu.rds'

library(Seurat)
library(mclust)

classifier <- c(
    "M1 Macrophage"              = 'cd68',
    "Helper T cells"             = 'CD4',
    "Cartilage tumor cell"       = 'sox9',
    "Proliferation Marker"       = 'ki67',
    "Tumor/Myeloid"              = 'lcp1',
    "MONOCYTE"                   = 'CD14',
    "M2 MACROPHAGE"              = 'cd206',
    "M2 MACROPHAGE"              = 'CD163',
    "HSC"                        = 'cd45',
    "Cytotoxic T Cells"          = 'cd8',
    "mDC/cDC - Dendritic cells"  = 'CD11c',
    "B CELLS"                    = 'cd19',
    "Endothelial"                = 'cd34',
    "Pericyte/Fibroblast"        = 'asma'
)

pilot <- readRDS(args[1])

dat <- GetAssayData(pilot, assay = 'RNA', layer = 'counts')
dat <- apply(dat, 1, function(marker){
    id <- order(marker) 
    dens <- densityMclust(marker[id], plot = FALSE, verbose = FALSE)
    cdf <- cdfMclust(dens, data=marker[id])
    cdf$y[match(marker, marker[id])]
}, simplify = FALSE)
dat <- do.call(rbind, dat)
colnames(dat) <- colnames(pilot)
pilot$GMM <- CreateAssayObject(data = dat)
dat <- dat>0.5
celltype <- dat[rownames(dat) %in% classifier, , drop=FALSE]
mode(celltype) <- 'character'
for(n in rownames(celltype)){
    celltype[n, ] <- ifelse(dat[n, ], names(classifier)[classifier==n], NA)
}
celltype <- celltype[match(classifier, rownames(celltype)), , drop=FALSE]
rownames(celltype) <- classifier
ct <- apply(celltype, 2, function(.e) {
    .w <- which(!is.na(.e))
    if(length(.w)>0){
        .e[.w[1]]
    }else{
        'unknown'
    }
})
stopifnot(identical(colnames(pilot), names(ct)))
write.csv(ct, file.path(dirname(args[1]), 'celltype.csv'))
pilot$celltype <- ct
saveRDS(pilot, sub('.rds$', '.withCellType.rds', args[1]))


