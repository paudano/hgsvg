library(getopt)
options <- matrix(c("inh", "i", 2, "character",
                    "pal", "p", 2, "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

t <- read.table(args$inh,header=T,comment.char="")

hom <- which(t$hap=="HOM")
h0 <- which(t$hap=="HAP0")
h1 <- which(t$hap=="HAP1")
maxCov <- 50

lc <- which(t$nAlt < 25 & t$nAlt.1 < 25)
altIdx <- intersect(lc, c(h0,h1))
homIdx <- intersect(lc,hom)

altTab <-table(t$nAlt[altIdx], t$nAlt.1[altIdx])
homTab <-table(t$nAlt[homIdx], t$nAlt.1[homIdx]) 



library(RColorBrewer)
library(lattice)

reds <- brewer.pal(9,args$pal)
par(mfrow=c(2,1))
levelAlt <- levelplot(altTab, col.regions=reds, at=seq(0,max(altTab), by=max(altTab)/length(reds)), xlab="Paternal coverage", ylab="Maternal coverage", main=sprintf("%s Het",args$sample), scales=list(x=list(at=c(seq(0,maxCov,by=5))),y=list(at=c(seq(0,maxCov,by=5)))))
levelHom <- levelplot(homTab, col.regions=reds, at=seq(0,max(homTab), by=max(homTab)/length(reds)), xlab="Paternal coverage", ylab="Maternal coverage", main=sprintf("%s Hom",args$sample), scales=list(x=list(at=c(seq(0,maxCov,by=5))),y=list(at=c(seq(0,maxCov,by=5)))))


pdf(sprintf("%s.parental_coverage.pdf",args$sample),height=4,width=8)
print(levelAlt,split=c(1,1,2,1),more=T)
print(levelHom,split=c(2,1,2,1),more=T)

dev.off()
