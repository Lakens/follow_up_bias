# Four conditions:
# 1 and 55: 2 groups, npilot = 10, small and large ES (and eta; 80% power)
# 25 and 79: 4 groups, npilot = 50, small and large ES (and eta; 80% power)
setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")
options(scipen=20) #disable scientific notation

cond <- 1

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,10)
n <- c(10000*(10:1))
for(i in 1:10){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES1 <- meanES
meanp1 <- meanp
nrNA1 <- nrNA

####

cond <- 55

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,10)
n <- c(10000*(10:1))
for(i in 1:10){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES55 <- meanES
meanp55 <- meanp
nrNA55 <- nrNA

####

cond <- 25

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,10)
n <- c(10000*(10:1))
for(i in 1:10){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES25 <- meanES
meanp25 <- meanp
nrNA25 <- nrNA

####

cond <- 79

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,10)
n <- c(10000*(10:1))
for(i in 1:10){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES79 <- meanES
meanp79 <- meanp
nrNA79 <- nrNA


##### General settings for plot

library("reshape2")
library("ggplot2")
library("scales")
cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cutoff<-seq(100,10,-10) #set cutoffs (! in reversed order, to match code above)

##### Plot p-mean
plotmeans<-data.frame(cutoff,meanES1, meanES55, meanES25, meanES79) #combine data into one dataframe
colnames(plotmeans) <- c("cutoff","k = 2, npilot = 10, ES = 0.0099","k = 2, npilot = 10, ES = 0.1379", "k = 4, npilot = 50, ES = 0.0099", "k = 4, npilot = 50, ES = 0.1379") #2 groups, npilot = 10, small and large ES (and eta; 80% power)
plotmeans <- melt(plotmeans, id="cutoff") #convert to long format with reshape package to easily plot multiple groups
colnames(plotmeans) <- c("cutoff", "condition","value") #change names in final dataframe to plot

tiff(file=paste("pmean.tiff",sep=""),width=3000,height=1500, units = "px", res = 300)
plot.new()
ggplot(data=plotmeans, aes(x = cutoff, y=value, group=condition, colour=condition)) +
  scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(size=1.5) + geom_point(size=3) + #plot lines and dots
  ylab("Mean p-value")  + xlab("cut-off (x 1000") + #set labels
  scale_x_continuous(breaks=pretty(plotmeans$cutoff, n = 10)) + #use scales package to set breaks
  guides(col = guide_legend(nrow = 2, byrow = TRUE)) + #split legend in two columns
  theme_bw(base_size=20) + theme(legend.position = "top") #increase font size, put legend on top
dev.off()


##### Plot NA
plotNA<-data.frame(cutoff,nrNA1, nrNA55, nrNA25, nrNA79) #combine data into one dataframe
colnames(plotNA) <- c("cutoff","k = 2, npilot = 10, ES = 0.0099","k = 2, npilot = 10, ES = 0.1379", "k = 4, npilot = 50, ES = 0.0099", "k = 4, npilot = 50, ES = 0.1379") #2 groups, npilot = 10, small and large ES (and eta; 80% power)
plotNA <- melt(plotNA, id="cutoff") #convert to long format with reshape package to easily plot multiple groups
colnames(plotNA) <- c("cutoff", "condition","value") #change names in final dataframe to plot

tiff(file=paste("NRexcludedvar.tiff",sep=""),width=3000,height=1500, units = "px", res = 300)
plot.new()
ggplot(data=plotNA, aes(x = cutoff, y=value, group=condition, colour=condition)) +
  scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(size=1.5) + geom_point(size=5) + #plot lines and dots
  ylab("nr excluded values")  + xlab("cut-off (x 1000)") + #set labels
  scale_x_continuous(breaks=pretty(plotNA$cutoff, n = 10)) + #use scales package to set breaks
  guides(col = guide_legend(nrow = 2, byrow = TRUE)) + #split legend in two columns
  theme_bw(base_size=20) + theme(legend.position = "top") #increase font size, put legend on top
dev.off()
