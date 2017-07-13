# Motivation of the number of replications
# See second part of the Supplementary Material

cond <- 1
# n contains the 1 million computed values of n(main)
# ES contains the 1 million computed effect sizes
n  <- readRDS(paste("n_",cond,".rds",sep=""))
ES <- readRDS(paste("ES_",cond,".rds",sep=""))

# Convert from 1 million values into 1000 blocks of 1000 values
nMAT <- matrix(n,nrow=1000,byrow=TRUE)
ESMAT <- matrix(ES, nrow=1000, byrow=TRUE)

# For each of the 1000 blocks, compute the estimated n and ES
fmean <- function(x){ mean(x,na.rm = TRUE) }
n_per1000 <- apply(nMAT, 1, fmean)
ES_per1000 <- apply(ESMAT, 1, fmean)

# CI for Effect size
round(mean(ES_per1000) - qt(p=.975,df=999) * sd(ES_per1000)/sqrt(1000),5)
round(mean(ES_per1000) + qt(p=.975,df=999) * sd(ES_per1000)/sqrt(1000),5)

# CI for n
mean(n_per1000) - qt(p=.975,df=999) * sd(n_per1000)/sqrt(1000)
mean(n_per1000) + qt(p=.975,df=999) * sd(n_per1000)/sqrt(1000)
