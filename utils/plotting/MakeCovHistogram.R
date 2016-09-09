args <- commandArgs(trailingOnly = TRUE)

tn = args[1]
on = args[2]
t <- read.table(tn)
pdf(on)
sd <- sd(t$V1)
m <- mean(t$V1)
md <- median(t$V1)
i <- which(t$V1 < m+3*sd & t$V1 > 4)
tl <- sprintf("%s mean %2.2f s.d. %2.2f median %d", tn, m, sd, md)
nbreaks=ceiling(m+3*sd)
par(xpd=T)
hist(t$V1[i], breaks=nbreaks,xlab="Coverage", ylab="Count", main=tl)
dev.off()

