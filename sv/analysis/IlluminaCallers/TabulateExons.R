#setwd("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/HG00514")
library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "op", "p",2, "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)
#args <- data.frame(tab="HG00514.merged_nonredundant_exon.DEL.bed", op="DEL", sample="HG00514")
t <- read.table(as.character(args$tab),header=T,comment.char="")

tab <- table(t$union, t$svClass)
library(RColorBrewer)


names <- rownames(tab)

methods <- c("PacBio", "Illumina", "BioNano", "PacBio,Illumina", "PacBio,BioNano", "Illumina,BioNano", "All")
for (i in 1:length(methods)) {
  if (length(which(names == methods[i])) == 0) {
    last <- length(tab[,"FRAMESHIFT"])
    tab <- rbind(tab, c(0,0,0,0))
    rownames(tab)[last+1] <- methods[i];
  }
}
p <- brewer.pal(5,"Set1")
pdf(sprintf("%s.%s.exon_class.pdf",args$sample, args$op),height=4,width=8)
barplot(t(tab[methods,]), col=p, main=sprintf("%s %s",args$sample, args$op))
legend("topright", legend=colnames(tab), pt.bg=p, pch=22)
dev.off()

martab <- rbind(tab, apply(tab,2,sum))
lastrow <- dim(martab)[1]
rownames(martab)[lastrow] <- "Total"

useRowNames <- T
if (args$op == "INS") {
  useRowNames <- F;
}

martab <- cbind(martab, apply(martab,1,sum))
lastcol <- dim(martab)[2]
colnames(martab)[lastcol] <- "Total"
tableName <- sprintf("%s.%s.exon_class.tsv",args$sample, args$op)
print("table name")
print(tableName)
order <- c("All","BioNano","Illumina","Illumina,BioNano","PacBio","PacBio,BioNano","PacBio,Illumina")
write.table(martab[order,], tableName, quote=F,row.names=useRowNames, sep="\t")
