library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "operation", "o", 2, "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)


args <- getopt(options)


outName <- paste(args$sample, "merged_nonredundant", args$operation, "bed.summary", sep=".")

t <- read.table(args$tab,header=T,comment.char="")

mergeTypes <- c("PacBio", "Illumina", "BioNano", "PacBio,Illumina", "PacBio,BioNano", "Illumina,BioNano", "All")
avgLen <- sapply(mergeTypes, function(i) mean(t$svLen[which(t$union == i)]))
n <- sapply(mergeTypes, function(i) length(which(t$union == i)))

#setwd("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/HG00733")
#t <- read.table("HG00733.merged_nonredundant.DEL.bed", header=T,comment.char="")

df <- data.frame(AverageLength=avgLen, Count=n, Op=rep(args$op, length(n)), MergeType=mergeTypes)

#df <- data.frame(AverageLength=avgLen, Count=n, Op=rep("DEL", length(n)), MergeType=mergeTypes)

#df[c("Op", "MergeType", "Count", "AverageLength")]

fdf <- format(df[c("Op", "MergeType", "Count", "AverageLength")], digits=5)
write.table(fdf, row.names=F, sep="\t", quote=F, file=outName)



