library(getopt)
options <- matrix(c("tab", "t", 2, "character",
                    "operation", "o", 2, "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

t <- read.table(args$tab,header=T,comment.char="")

pass <- which(t$orth_filter ==  "PASS")
fail <- which(t$orth_filter ==  "FAIL")


short <- 1000
long <- 20000
minLength <- 50
shorti <- which(t$svLen >= 50 & t$svLen <= short )
longi  <- which(t$svLen >= short & t$svLen <= long )

shortPass <- intersect(shorti, pass)
shortFail <- intersect(shorti, fail)
hShortPass <- hist(t$svLen[shortPass], breaks=seq(50,short,by=50),plot=F)
hShortFail <- hist(t$svLen[shortFail], breaks=seq(50,short,by=50),plot=F)





library(RColorBrewer)
outFile <- paste(args$sample,args$operation,"pass","pdf", sep=".")

print(outFile)

pdf(outFile, width=12,height=6)
par(mfrow=c(1,2))
ymax <- max(c(hShortPass$counts, hShortFail$counts))
p <- brewer.pal(3,"Set1")

plot(c(), xlim=c(0,short+50), ylim=c(0,ymax*1.1), xlab="SV Length", ylab="Count",main=paste(args$sample, args$operation, "SVs", "\n", "50-1kbp", sep=" " ))
sapply(seq(1,length(hShortPass$counts)), function(i) rect(hShortPass$mids-25, 0, hShortPass$mids-10, hShortPass$counts, col=p[2]))
sapply(seq(1,length(hShortFail$counts)), function(i) rect(hShortFail$mids-5, 0, hShortFail$mids+10, hShortFail$counts, col="grey"))
legend("topright", legend=c("Pass", "Fail"), pch=22,pt.bg=c(p[2], "grey"),pt.cex=1.5)



longPass <- intersect(longi, pass)
longFail <- intersect(longi, fail)
hLongPass <- hist(t$svLen[longPass], breaks=seq(short,long,by=1000), plot=F)
hLongFail <- hist(t$svLen[longFail], breaks=seq(short,long,by=1000),plot=F)

ymax <- max(c(hLongPass$counts, hLongFail$counts))


p <- brewer.pal(3,"Set1")

plot(c(), xlim=c(0,long+50), ylim=c(0,ymax*1.1), xlab="SV Length", ylab="Count",main=paste(args$sample, args$operation, "SVs", "\n", "1kbp-20kbp", sep=" " ))
sapply(seq(1,length(hLongPass$counts)), function(i) rect(hLongPass$mids-500, 0, hLongPass$mids-100, hLongPass$counts, col=p[2]))
sapply(seq(1,length(hLongFail$counts)), function(i) rect(hLongFail$mids-50, 0, hLongFail$mids+350, hLongFail$counts, col="grey"))
legend("topright", legend=c("Pass", "Fail"), pch=22,pt.bg=c(p[2], "grey"),pt.cex=1.5)

dev.off()






outFile <- paste(args$sample,args$operation,"line","pdf", sep=".")
fp <- brewer.pal(4,"Paired")
print(outFile)

pdf(outFile, width=12,height=6)
par(mfrow=c(1,2))
shortPass <- intersect(shorti, pass)
shortFail <- intersect(shorti, fail)
par(mar = c(5,5,2,5))
hShortPass <- hist(t$svLen[shortPass], breaks=seq(50,short,by=50),col=fp[1], main=paste(args$sample, args$operation, "SVs", "\n", "50-1kbp", sep=" " ), xlab="SV length")
hShortFail <- hist(t$svLen[shortFail], breaks=seq(50,short,by=50),plot=F)
par(new=T)
plot(hShortPass$mids, hShortPass$counts/(hShortPass$counts+hShortFail$counts), axes=F, xlab=NA, ylab=NA, pch=21, bg=fp[2], ylim=c(0,1))
axis(side = 4)

mtext(side=4,line=2,"Fraction of validated sites")



longPass <- intersect(longi, pass)
longFail <- intersect(longi, fail)
par(mar = c(5,5,2,5))
hLongPass <- hist(t$svLen[longPass], breaks=seq(short,long,by=1000), plot=T, col=fp[1],main=paste(args$sample, args$operation, "SVs", "\n", "1kbp-20kbp", sep=" " ), xlab="SV length")
hLongFail <- hist(t$svLen[longFail], breaks=seq(short,long,by=1000),plot=F)

par(new=T)
plot(hLongPass$mids, hLongPass$counts/(hLongPass$counts+hLongFail$counts), axes=F, xlab=NA, ylab=NA, pch=21, bg=fp[2], ylim=c(0,1))
axis(side = 4)

dev.off()

