#!/usr/bin/env Rscript
# This script allows to calculate the dynamic range.
# The standard approach assessing dynamic range is to calculate a signal-to-background (SNR) ratio by dividing the average of the top 20 brightest cells by the average intensity of the weakest 10% of cells.
# An SNR of 10 or more supports reliable image analysis.
# The recommend SNR range is (10, Inf], typically > 100.
# Lower than 3 indicates pool-performing antibodies.
# 
# Script author: Jianhong Ou
# Script should work after running ExportCellDetectionMeasurement.groovy
# This source code is licensed under the MIT license
# How to run in command line:
# Rscript dynamicRanges.R path/to/measurementsExport/folder
# How to run in Rstudio/R:
# uncomment line 18 and change the folder to the measurements folder
# Output: Inf indicates zero background. NaN indicates zero/zero
#         Column names is the measurement method for each cell
args = commandArgs(trailingOnly=TRUE)
# args <- c(path='~/Downloads/qupath/run_20230926/measurementsExport13635325425080860669')
fs <- dir(args[1], "qptiff.tsv", full.names=TRUE)
useValue <- c('Mean', 'Median', 'Min', 'Max')

## meta data information filter
metaInfoStartWith <- '^(Nucleus|Cell)'

for(f in fs){
    message('Working on:', f)
    pf <- dirname(f)
    prefix <- sub(".qptiff.tsv", '', basename(f))
    d <- read.delim(f)
    d <- d[, !grepl(metaInfoStartWith, colnames(d))]
    cn <- colnames(d)
    genes <- sub("..Cell..Median", '', cn[grepl('Cell..Median', cn)])
    dx <- lapply(useValue, function(.ele){
        dat <- d[, cn[grepl(paste0('Cell..', .ele), cn)]]
        colnames(dat) <- sub(paste0('..Cell..', .ele), '', colnames(dat))
        dat <- t(dat)
        colnames(dat) <- rownames(meta)
        dat
    })
    names(dx) <- useValue
    top20 <- seq.int(20)
    t20 <- lapply(dx, function(.ele){
        n <- ncol(.ele)
        apply(.ele, 1, function(x){
            x <- sort(x, decreasing = TRUE)
            ## top 20 vs bottom 10
            mean(x[top20], na.rm=TRUE)
        }, simplify = TRUE)
    })
    t20 <- do.call(cbind, t20)
    write.csv(t20, file.path(pf, paste0(prefix, '.top20cells.csv')))
    dr <- lapply(dx, function(.ele){
        n <- ncol(.ele)
        not_bottom10 <- seq.int(floor(n*.9))
        apply(.ele, 1, function(x){
            x <- sort(x, decreasing = TRUE)
            ## top 20 vs bottom 10
            mean(x[top20], na.rm=TRUE)/mean(x[-not_bottom10], na.rm=TRUE)
        }, simplify = TRUE)
    })
    dr <- do.call(cbind, dr)
    write.csv(dr, file.path(pf, paste0(prefix, '.dynamicRange.top20cells.vs.bottom10percentCells.csv')))
}