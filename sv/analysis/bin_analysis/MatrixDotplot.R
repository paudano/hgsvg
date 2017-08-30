library(getopt)
options <- matrix(c("dots", "d", 2, "character",
                    "region", "r", 2, "character",
                    "tmpdir", "t", 2, "character",
                    "output", "o", 2, "character"),
                  byrow=T, ncol=4)

args <- getopt(options)
sources <- strsplit(args$dots, " ")[[1]]
ncol <- length(sources)
pdf(args$output)
print(sources)
print(args$tmpdir)

layoutMat <- layout(matrix(1:(ncol^2), ncol, ncol, byrow=T), widths=c(1,rep(4,ncol-1)), heights=c(1,rep(4,ncol-1)))

par(mar=c(1,0,0,0))      
plot.new()
pn <- 1
for (i in seq(2,ncol)) {
  par(mar=c(1,0,0,0))      
  plot.new()
  mtext(sources[i],side=1,cex=0.75);
  pn <- pn+1
}

for (i in seq(1,(ncol-1))) {
  # add marginal title
  par(mar=c(0,2,0,0))      
  plot.new()
  
  pn <- pn+1
  for (j in seq(2,ncol)) {
   if (j <= i) {
     par(mar=c(1,0,0,0))
     plot.new();
     if (i == ncol-1 & j == 3) {
       if (is.null(args$region) == F) {
         mtext(args$region,side=1,cex=1);
       }
     }
   } else {
      dotsFileName <- paste(args$tmpdir, paste(sources[i], sources[j], "dots", sep="."), sep="/")
      t <- read.table(dotsFileName)
      xRange <- range(t$V1)
      xRange[1] <- max(0,xRange[1])
      yRange <- range(t$V2)
      par(mar=c(1,1,1,1))
      
      plot(c(),xlim=xRange,ylim=yRange,xlab="",ylab="",  main="", cex.axis=0.75)
      if (j == i+1) {
        mtext(sources[i], side=2,line=2,cex=0.75);
      }
      
      segments(t$V1, t$V2, t$V1+t$V3, t$V2+t$V3, col="black", lwd=0.5);
      pn <- pn+1;
    }
   
#    if (i == 1) {
#      mtext(sources[j],side=3, line=2, cex=0.75)
#    }
    
  }
}


dev.off()
