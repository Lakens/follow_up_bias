# Detailed look at main example, K = 2, 

setwd("C:\\Users\\Daniel\\surfdrive\\Word\\ESpaper Casper\\Power Casper")

# Check which condition numbers simulated which data:
K      <- c(2, 3, 4)        # number of groups in a one-way anova
npilot <- c(10, 25, 50)     # sample size per group in the pilot
ES     <- c(.0099,.0588,.1379) # population ES
alpha  <- 0.05              # nominal level of significance
beta   <- c(0.2, 0.1)       # 1 - power
ESmeasure <- c("eta","epsilon","omega")
R      <- 1000000              # number of replications per condition. Set this to 10^6 at the end
conditions <- expand.grid(ESmeasure,K, npilot, ES, alpha, 1-beta)
conditions #see the list of 162 conditions

nstar <- seq(250,10)  # note: the values *must* be listed in decreasing order
conditions <- c(10, 37, 64) #N = 25
#conditions <- c(1, 28, 55) #N = 10
#conditions <- c(19, 46, 73) # N = 50
powermatrix <- matrix(NA, nrow=length(conditions),ncol=length(nstar))
colnames(powermatrix) <- nstar
rownames(powermatrix) <- conditions

for(x in 1:length(conditions)){ # x = condition
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

powermat<-as.matrix(powermatrix)
powermat<-t(powermat)
colnames(powermat) <- c("ES = 0.0099","ES = 0.0588","ES = 0.1379") 

#Get the max power when n* is 100,000
max_obs_power <- numeric(length(conditions))
for(x in 1:length(conditions)){ # x = condition
  cond <- conditions[x]
  pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))
  max_obs_power[x] <- sum(pvalueMAT<0.05, na.rm=TRUE)/1000000
}

library("reshape2")
powermat <- melt(powermat, id="ES") #convert to long format with reshape package to easily plot multiple groups
colnames(powermat) <- c("N","ES","Power") 
#cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette<-c("#5e3c99", "#e66101", "#fdb863")


library("ggplot2")
tiff(file=paste("Fig7_",toString(conditions),".tiff",sep=""),width=2700,height=1700, units = "px", res = 300)
ggplot(data=powermat, aes(x = N, y=Power, group=ES)) +
  geom_hline(yintercept = max_obs_power[1], size=1, alpha=0.4, colour = "black" ) +
  geom_hline(yintercept = max_obs_power[2], size=1, alpha=0.4, colour = "black", linetype = 2) +
  geom_hline(yintercept = max_obs_power[3], size=1, alpha=0.4, colour = "black", linetype = 3) +
  scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(aes(linetype=ES),size=1.5) + #geom_point(size=5) + #plot lines and dots
  scale_linetype_manual(values=c(1,2,3)) + 
  ylab("Power in follow-up study")  + xlab("Maximum sample size per condition") + #set labels
  guides(col = guide_legend(nrow = 2, byrow = TRUE)) + #split legend in two columns
  theme_bw(base_size=20) + expand_limits(y=c(0,1)) + theme(legend.title=element_blank(), legend.position = "top", legend.key.width=unit(2,"cm")) #increase font size, put legend on top
dev.off()


