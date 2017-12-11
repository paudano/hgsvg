library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "op", "p",2, "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

setwd("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/HG00514")

args <- data.frame(tab="tmp.txt", sample="HG00514", op="DEL")
t <- read.table(as.character(args$tab),header=T,comment.char="")
head(t)
pbi <- which(as.character(t$union) == "PacBio" |
             as.character(t$union) == "PacBio,Illumina" |
             as.character(t$union) == "PacBio,BioNano" |
             as.character(t$union) == "All")

ili <- which(as.character(t$union) == "Illumina" |
             as.character(t$union) == "PacBio,Illumina" |
             as.character(t$union) == "Illumina,BioNano" |
             as.character(t$union) == "All")


bni <- which(as.character(t$union) == "BioNano" |
             as.character(t$union) == "PacBio,BioNano" |
             as.character(t$union) == "Illumina,BioNano" |
             as.character(t$union) == "All")

head(t)

npb <- length(which(t$pc_frac[pbi] > 0))
nil <- length(which(t$pc_frac[pbi] > 0))
nbn <- length(which(t$pc_frac[bni] > 0))

