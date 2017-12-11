library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "op", "p",2, "character",
                    "method", "m",2, "character",                    
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

t <- read.table(args$tab,header=T,comment.char="")


library(RColorBrewer)
pal <- brewer.pal(7,"Set1")

pbc <- which(t$callset=="pacbio")
bnc <- which(t$callset=="bionano")
ilc <- which(t$callset=="Illumina")

pbi <- which(t$in_pb == "TRUE")
bni <- which(t$in_bn == "TRUE")
inti <- which(t$in_int == "TRUE")
if (args$op == "DEL") {
  opTitle="Deletion";
} else {
  opTitle="Insertion";
}

pbOnlyi <- intersect(setdiff(pbi, union(bni,inti)), pbc)
bnOnlyi <- intersect(setdiff(bni, union(pbi,inti)), bnc)
intOnlyi <- intersect(setdiff(inti, union(pbi,bni)), ilc)

pbAndBn <- intersect(setdiff(intersect(pbi,bni), inti), pbc)
pbAndIl <- intersect(setdiff(intersect(pbi,inti), bni), pbc)
bnAndIl <- intersect(setdiff(intersect(bni,inti), pbi), bnc)

allThree <- intersect(intersect(intersect(bni,inti),pbi),pbc)

shortMax <- 1000
short <- which(t$svLen <= shortMax & t$svLen >= 50)
longMax <- 20000
long <-  which(t$svLen >= shortMax & t$svLen < longMax)

idxList <- list(pbOnlyi,bnOnlyi, intOnlyi, pbAndBn, pbAndIl, bnAndIl, allThree)

shortIdxList <- sapply(idxList, function(x) intersect(short, x))

shortBreaks <- seq(50,shortMax,by=50)
print("short hist")
shortHist <- lapply(shortIdxList, function(x) hist(t$svLen[x],breaks=shortBreaks,plot=F,right=F))

filename <- paste("IntegratedBarchart", args$sample, args$op, "pdf",sep=".")
print("writing to ")
print(filename)
pdf(filename,width=8,height=6)
par(mfrow=c(1,2))
colSum <- sapply(seq(1,length(shortBreaks)-1), function(i) sum(sapply(seq(1,length(shortHist)), function(hi) shortHist[[hi]]$counts[i])))
yStart <- rep(0,length(shortHist[[1]]$counts))

plot(c(), xlim=c(0,shortMax),ylim=c(0,max(colSum)), ylab="Count", xlab="Structural variant length (bp)", main=paste(args$method, "\n", args$sample, opTitle, " < 1000 bp", sep=" "))
for (hidx in seq(1,length(shortHist))) {
  lapply(seq(1,length(shortHist[[hidx]]$counts)), function(i) rect(shortHist[[hidx]]$mids[i]-25, yStart[i], shortHist[[hidx]]$mids[i]+25,yStart[i]+shortHist[[hidx]]$counts[i], col=pal[hidx]));
  yStart = yStart + shortHist[[hidx]]$counts;
  print(yStart)
}

methods <- c("PacBio", "BioNano", "Illumina", "PacBio,BioNano", "PacBio,Illumina", "BioNano,Illumina", "All")
legend("topright", legend=methods, pch=22,pt.bg=pal)



longIdxList <- sapply(idxList, function(x) intersect(long, x))

longBreaks <- seq(1000,longMax,by=1000)
longHist <- lapply(longIdxList, function(x) hist(t$svLen[x],breaks=longBreaks,plot=F,right=F))
colSum <- sapply(seq(1,length(longBreaks)-1), function(i) sum(sapply(seq(1,length(longHist)), function(hi) longHist[[hi]]$counts[i])))

yStart <- rep(0,length(longHist[[1]]$counts))

plot(c(), xlim=c(0,longMax),ylim=c(0,max(colSum)), ylab="Count", xlab="Structural variant length (bp)",main=paste(args$method, "\n", args$sample, opTitle, "\n >= 1000, < 20,000 bp", sep=" "))
for (hidx in seq(1,length(longHist))) {
  print(length(longIdxList[[hidx]]))
  lapply(seq(1,length(longHist[[hidx]]$counts)), function(i) rect(longHist[[hidx]]$mids[i]-500, yStart[i], longHist[[hidx]]$mids[i]+500,yStart[i]+longHist[[hidx]]$counts[i], col=pal[hidx]));
  yStart = yStart + longHist[[hidx]]$counts;
  print(yStart)
}

methods <- c("PacBio", "BioNano", "Illumina", "PacBio,BioNano", "PacBio,Illumina", "BioNano,Illumina", "All")
legend("topright", legend=methods, pch=22,pt.bg=pal)

dev.off()

