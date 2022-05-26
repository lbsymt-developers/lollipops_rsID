library(trackViewer)

load("data/cred_promoters_GRanges.rda")
hic <- readRDS("data/snps_mapping_HiC_granges.rds")
hic <- data.frame(hic)
hic <- hic[!duplicated(hic$rsid), ]

hic$rsid
credpromoter$rsid
