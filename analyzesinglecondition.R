setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")
options(scipen=20) #disable scientific notation

cond <- 37

nMAT  <- readRDS(paste("n_",cond,"N100max.rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,"N100max.rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,"N100max.rds",sep=""))

hist(ESMAT, breaks=400, xlim=c(0,0.1))
meanES <- mean(ESMAT,na.rm=TRUE)
meanES

obspower2 <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)

sum(pvalueMAT < .05, na.rm=TRUE)/(1000000-sum(is.na(pvalueMAT)))*100

cond <- 37

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

hist(ESMAT, breaks=400, xlim=c(0,0.1))
meanES <- mean(ESMAT,na.rm=TRUE)
meanES

sum(pvalueMAT < .05, na.rm=TRUE)/(1000000-sum(is.na(pvalueMAT)))*100
