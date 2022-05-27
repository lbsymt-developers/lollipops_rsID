library(trackViewer)

load("data/cred_promoters_GRanges.rda")
hic <- readRDS("data/snps_mapping_HiC_granges.rds")
hic <- data.frame(hic)
hic <- hic[!duplicated(hic$rsid), ]

credpromoter <- data.frame(credpromoter)
credpromoter <- credpromoter[!duplicated(credpromoter$start),]
load("data/promoter_ranges.rda")
promoters <- data.frame(promoterranges)


### LOLLIPOPS DE ABCA7

sample.rs4147929 <- GRanges("chr19", IRanges(start = c(1063444),
                                      end = c(1063444),
                                      width=1, names="rs4147929"))
sample.rs4147929$color <- "red"
sample.rs4147929$border <- "black"
sample.rs4147929$alpha <- 0.8
sample.rs4147929$label.parameter.rot <- 45

sample.rs3752241 <- GRanges("chr19", IRanges(start = c(1053524),
                                             end = c(1053524),
                                             width=1, names="rs3752241"))
sample.rs3752241$color <- "red"
sample.rs3752241$border <- "black"
sample.rs3752241$alpha <- 0.8
sample.rs3752241$label.parameter.rot <- 45

features <- GRanges("chr19", IRanges(c(1063051, 1064051),
                                    width=c(1000, 725),
                                    names=c("Promoter_region",
                                            "ENST00000531478")))
features$fill <- c("#FF8833", "#51C6E6")
# features$height <- c(0.02, 0.03)

features_2 <- GRanges("chr19", IRanges(c(1053075, 1054076),
                                     width=c(1000, 571),
                                     names=c("Promoter_region",
                                             "ENST00000530092")))
features_2$fill <- c("#FF8833", "#51C6E6")

gene <- GRanges("chr19", IRanges(c(1039997),
                                 width=c(25576),
                                 names=c("ABCA7")))
gene$fill <- "#51C6E6"

snps <- GRanges("chr19", IRanges(start = c(1053524, 1063444),
                                 end = c(1053524, 1063444),
                                 width=1,
                                 names=c("rs3752241", "rs4147929")))
snps$color <- "red"
snps$border <- "black"
snps$alpha <- 0.8
snps$label.parameter.rot <- 45

# FACTORES TRANSCIRPCIONALES EN LA REGION

peak_scores <- c(21.4, 17.4, 5.13, 9.12, 11.6, 3.27, 14.1)
ft <- GRanges("chr19", IRanges(start = c(1039887, 1038804,
                                         1039389, 1039399,
                                         1039786, 1040018,
                                         1039961),
                               end = c(1040214, 1040306,
                                       1039503, 1040259,
                                       1040317, 1040160,
                                       1040190),
                               names = c("NRF1", "POLR2A",
                                         "REST", "SIN3A",
                                         "SPI1", "TAF1",
                                         "USF1")
                               ))
cols <- palette("Classic Tableau")
ft$fill <- cols[1:7]

tiff("images/lollipop_ABCA7.tiff", height = 20, width = 20, units='cm',
     compression = "lzw", res = 300)
lolliplot(list(A = sample.rs4147929, B = sample.rs3752241,
               C = snps),
          list(x = features, y = features_2, z = gene))
dev.off()

