# New study on effect of lowering n^*
# Two conditions: {K=2, npilot =10} (most difficult) and {K=4, npilot=50} (easiest)
# Population ES either S, M or L
# Only focus on eta^2 (because other 2 are similar; and we want to keep things overseeable)
# In total: 6 conditions
setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")

# Focus on following n^*:
nstar <- c(100000, 10000, 1000, 200, 100,50)  # note: the values *must* be listed in decreasing order
conditions <- c(1:162)   #
powermatrix <- matrix(NA, nrow=length(conditions),ncol=length(nstar))
colnames(powermatrix) <- nstar
rownames(powermatrix) <- conditions

for(x in 1:162){ # x = condition
  cond <- conditions[x]
  nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
  ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
  pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))
  powerr <- meanES <- meanp <- nrNA <- rep(NA,5)
  for(i in 1:length(nstar)){
    # Same as in 'loweringN.R'
    kill <- (nMAT> nstar[i])
    pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
    meanES[i] <- mean(ESMAT,na.rm=T)
    meanp[i] <- mean(pvalueMAT,na.rm=T)
    nrNA[i] <- sum(is.na(nMAT))  
    # Power: proportion of 'significant' results out of all results
    # Count those with p<.05 relative to all (and all = 'p < 2' :-) )
    powerr[i] <- sum(pvalueMAT<.05, na.rm=T)/(sum(pvalueMAT < 2, na.rm=T))
  }
  powermatrix[x,] <- powerr
}

write.table(powermatrix, file="LowerMaxNTableS1.csv")