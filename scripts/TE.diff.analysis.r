#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(xtail))
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-m", "--mrna", required=TRUE, help="input mRNA counts")
parser$add_argument("-r", "--rpf", required=TRUE, help="input RPF counts")
parser$add_argument("-o", "--output", required=TRUE, help="output TE diff table")
parser$add_argument("-p", "--case", required=TRUE, help="case sample ids")
parser$add_argument("-n", "--control", required=TRUE, help="control sample ids")
args <- parser$parse_args()

message("load RNA-seq counts ...")
RNA.seq.counts <- read.table(args$mrna,sep="\t",row.names = 1,header = T, check.names = F)
message("load ribo-seq counts ...")
ribo.seq.counts <- read.table(args$rpf,sep="\t",row.names = 1,header = T, check.names = F)

case.ids <- unlist(strsplit(args$case,","))
control.ids <- unlist(strsplit(args$control,","))

message(paste0("case samples: ", args$case))
message(paste0("control samples: ", args$control))

sample.ids <- c(case.ids,control.ids)
conditions <- c(rep("case",length(case.ids)),rep("control",length(control.ids)))
mrna <- RNA.seq.counts[,sample.ids] + 1
rpf <- ribo.seq.counts[,sample.ids] + 1

message("start xtail analysis ...")
res <- xtail(mrna, rpf, conditions, baseLevel="control", minMeanCount=1, bins=5000)

message(paste0("saving results to ", args$output))
write.table(res$resultsTable, args$output ,sep="\t", quote=F)
message("all done.")
