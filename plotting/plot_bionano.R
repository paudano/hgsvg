library(getopt)
options <- matrix(c("pbsv", "p", 2, "character",
                    "bndel", "d", 2, "character",
                    "bnins", "i", 2, "character",                    
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)


t <- read.table(args$pbsv,header=TRUE, comment.char="")

                                        #
# First plot the confirmed calls
ratioLimit <- 0.25

del <- which(t$bnOverlap < ratioLimit  & t$svType=="deletion")
ins <- which(t$bnOverlap < ratioLimit  & t$svType=="insertion")
indels <- c(del,ins)
maxx <- max(c(length(del),length(ins)))
maxy <- max(c(t$svLen[del], t$svLen[ins]))

ptcol=rep("black", length(c(del,ins)))
confdel <- which(t$svType[indels] == "deletion")
ptcol[confdel] = "red"
options(scipen=20)
sn <- paste("PacBioConfirmedCalls.",args$sample,".pdf",sep="")

pdf(sn)
print(args)
m <- sprintf("Confirmed %s", args$sample)
plot(sort(t$svLen[del]),pch=16, col="red", xlim=c(0,maxx), ylim=c(100,maxy), ylab="SV Length", log="y", main=m)
points(sort(t$svLen[ins]),pch=16,col="black")
legend("topleft", legend=c("Insertion", "Deletion"),pch=16,col=c("black", "red"))
dev.off()

# Investigate confirmation by length
bins <- seq(1000, 20000,by=1000)
delCount <- sapply(bins, function(i) length(which(t$svType=="deletion" & t$svLen > i)))
confDelCount <- sapply(bins, function(i) length(which(t$svType=="deletion" & t$svLen > i & t$bnOverlap < ratioLimit)))

insCount <- sapply(bins, function(i) length(which(t$svType=="insertion" & t$svLen > i)))
confInsCount <- sapply(bins, function(i) length(which(t$svType=="insertion" & t$svLen > i & t$bnOverlap < ratioLimit)))

write.table(matrix(c(bins, delCount, confDelCount, insCount, confInsCount), ncol=5),row.names=F, col.names=c("SV Size", "Del count", "Conf. del", "Ins. count", "Conf. ins. count"), sep=",", file=sprintf("BioNanoTable.%s.csv",args$sample))

library(RColorBrewer)
bl <- brewer.pal(7,"Blues")
rd <- brewer.pal(7,"Reds")
gy <- brewer.pal(7,"Greys")

pdf(sprintf("FractionValidatedSVs.%s.pdf",args$sample))
par(mar=c(5,4,4,4)+0.1)
plot(bins, confDelCount / delCount, pch=21, bg=rd[4],cex=2, ylab="Fraction of validated BioNanoGenomics calls", xlab="Total calls >= size", main=args$sample)
par(new=TRUE)
plot(bins, delCount, axes=F,ylab="",xlab="", col="black", bg=rd[7], pch=22,cex=2, ylim=c(0,max(c(delCount,insCount))))
par(new=TRUE)
plot(bins, confInsCount / insCount, pch=21, bg=gy[4], xlab="",ylab="", axes=F,cex=2)
par(new=TRUE)
plot(bins, insCount, axes=F,ylab="",xlab="", col="black", bg=gy[7], pch=22,cex=2, ylim=c(0,max(c(delCount,insCount))))

axis(4, seq(0,max(c(insCount, delCount)), by=100))
mtext("Number of SVs", side=4, line=2)
legend("bottomleft", c("Validating del.", "Total del.", "Valid ins.", "Total ins."), pch=c(21,22,21,22), pt.bg=c(rd[4], rd[7], gy[4], gy[7]))
dev.off()


bnDel <- read.table(args$bndel, comment.char="", header=T)
bnIns <- read.table(args$bnins, comment.char="", header=T)




ni <- which(bnDel$svLen_2 <= 1)
bnSvLen <- c(bnDel$svLen, bnIns$svLen)
pbSvLen <- c(bnDel$svLen_2, bnIns$svLen_2)

pdf(sprintf("BioNano.PacBio.Agreement.%s.pdf",args$sample))
i <- which(bnDel$svLen_2 > 1)
plot(bnDel$svLen[i], bnDel$svLen_2[i], log="xy", xlim=c(100,50000),ylim=c(100,50000),pch=16,cex=0.5,col="#FF000044", xlab="BioNanoGenomic SV Length", ylab="PacBio SV Length", main=sprintf("BioNanoSVVersusPacBio.%s.pdf",args$sample))
i <- which(bnIns$svLen_2 > 1)
points(bnIns$svLen[i], bnIns$svLen_2[i], pch=16,cex=0.5, col="#00000044")
legend("topleft", legend=c("Insertion", "Deletion"),pch=16,col=c("black", "red"))
dev.off()
