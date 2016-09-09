args <- commandArgs(trailingOnly = TRUE)

tn = args[1]
t <- read.table(tn)
sd <- sd(t$V1)
m <- mean(t$V1)
md <- median(t$V1)
i <- which(t$V1 < m+3*sd & t$V1 > 4)
tl <- sprintf("%s\t%2.2f\t%2.2f\t%d\n", tn, m, sd, md)
cat(tl)

