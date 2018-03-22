library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "operation", "o", 2, "character",
                    "ncaller", "n", 2, "integer",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)
#setwd("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/HG00514")
#args <- data.frame(tab="integrated.DEL.bed.passfail",operation="DEL", sample="HG00514")

if (is.null(args$ncaller) == TRUE) {
  args$ncaller=12
}

pf <- read.table(as.character(args$tab),header=T,comment.char="")
library(RColorBrewer)
p <- brewer.pal(3,"Set1")

ncallers <- sapply(pf$ALGORITHM, function(i) length(strsplit(as.character(i), split=",")[[1]]))
outName <- paste(as.character(args$sample), args$operation, "method_count.pdf", sep=".")

pft <- table(ncallers, pf$orth_filter)
pft
sd <- setdiff(seq(1,args$ncaller), as.numeric(rownames(pft)))
f <- as.numeric(pft[,"PASS"])/(as.numeric(pft[,"PASS"])+as.numeric(pft[,"FAIL"]))
f
pft <- cbind(pft, f)

pft
colnames(pft)[3] <- "Fraction"
print("after adding fraction")
pft
for (i in sd) {
  pft <- rbind(pft, c(0, 0, 0));
  idx <- dim(pft)[1];
  rownames(pft)[idx] <- i;
}
print("before plot")
pft
pdf(outName)
t(pft[,c("FAIL","PASS")])
barplot(t(pft[,c("FAIL","PASS")]),
        col=c("grey", p[2]),
        xlab="Number of algorithms",
        ylab="Number of calls",
        beside=T, main=paste(args$sample, args$operation, sep=" "))
legend("topright", legend=c("No evidence", "Pass"), pch=22, pt.bg=c("grey", p[2]))

pft <- cbind(pft, rep(args$operation, dim(pft)[1]))
colnames(pft)[4] <- "Op"

dev.off()
tableOutName <- paste(args$sample, args$operation, "method_count.txt", sep=".")
pft[is.nan(pft)] <- 0

write.table(pft,sep="\t", file=tableOutName, quote=F)
