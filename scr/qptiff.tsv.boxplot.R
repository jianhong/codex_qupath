#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)
library(reshape2)
fs <- dir(args[1], "*.qptiff.tsv", full.names = TRUE)
ylim <- c(0, as.integer(args[2]))
for(f in fs){
    x <- read.delim(f)
    m <- x[, grepl('Mean', colnames(x))]
    colnames(m) <- sub("..Mean", "", colnames(m))
    head(colnames(m))
    gp <- sub("^.*\\.\\.(Nucleus|Cytoplasm|Membrane|Cell)", "\\1", colnames(m))
    table(gp)
    m <- split(as.data.frame(t(m)), gp)
    m <- lapply(m, t)
    m1 <- lapply(m, melt)
    mapply(m1, names(m1), FUN=function(.ele, .n){
        .ele$Var2 <- sub(paste0("..", .n), "", .ele$Var2)
        p <- ggplot(.ele, aes(x=Var2, y=value, fill=Var2)) + geom_boxplot() + 
            ylim(ylim) +
            xlab("") + ylab(.n) + theme_minimal()
        ggsave(paste0(sub("qptiff.tsv", "", f), .n, ".pdf"), plot=p, width=24)
    })
}


