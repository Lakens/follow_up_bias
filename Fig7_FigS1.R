# Four conditions:
# 1 and 55: 2 groups, npilot = 10, small and large ES (and eta; 80% power)
# 25 and 79: 4 groups, npilot = 50, small and large ES (and eta; 80% power)
# 28 and 53: 2 groups, npilot = 10, small and large ES (and eta; 80% power)

setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")
options(scipen=20) #disable scientific notation

cond <- 1

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

n <- c(10000*(10:1), 1000*(9:1), 100*(9:1))
meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
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

cond <- 28

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES28 <- meanES
meanp28 <- meanp
nrNA28 <- nrNA

####

cond <- 53

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
  kill <- (nMAT> n[i])
  pvalueMAT[kill] <- ESMAT[kill] <- nMAT[kill] <- NA
  meanES[i] <- mean(ESMAT,na.rm=T)
  meanp[i] <- mean(pvalueMAT,na.rm=T)
  nrNA[i] <- sum(is.na(nMAT))  
}

meanES53 <- meanES
meanp53 <- meanp
nrNA53 <- nrNA

####

cond <- 55

nMAT  <- readRDS(paste("n_",cond,".rds",sep=""))
ESMAT <- readRDS(paste("ES_",cond,".rds",sep=""))
pvalueMAT <- readRDS(paste("p_",cond,".rds",sep=""))

meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
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

meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
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

meanES <- meanp <- nrNA <- rep(NA,length(n))
for(i in 1:length(n)){
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

cutoff<-n #set cutoffs (! in reversed order, to match code above)

##### Plot ES-mean
plotmeans<-data.frame(cutoff,meanES1, meanES28, meanES55, meanES25, meanES53, meanES79) #combine data into one dataframe
colnames(plotmeans) <- c("cutoff","k = 2, npilot = 10, ES = 0.0099","k = 2, npilot = 10, ES = 0.0588","k = 2, npilot = 10, ES = 0.1379","k = 4, npilot = 50, ES = 0.0099", "k = 4, npilot = 50, ES = 0.0588", "k = 4, npilot = 50, ES = 0.1379") #2 groups, npilot = 10, small and large ES (and eta; 80% power)
plotmeans <- melt(plotmeans, id="cutoff") #convert to long format with reshape package to easily plot multiple groups
colnames(plotmeans) <- c("cutoff", "condition","value") #change names in final dataframe to plot

tiff(file=paste("Fig7.tiff",sep=""),width=3800,height=2000, units = "px", res = 300)
plot.new()
ggplot(data=plotmeans, aes(x = cutoff, y=value, group=condition, colour=condition)) +
  scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(size=1.5) + # geom_point(size=5) + #plot lines and dots
  ylab("Mean ES-value")  + xlab("cut-off") + #set labels
  scale_x_log10(breaks=c(200,500,1000,2000,5000,10000,20000,50000,100000)) +
  # scale_x_continuous(breaks=pretty(plotmeans$cutoff, n = 10)) + #use scales package to set breaks
  guides(col = guide_legend(nrow = 2, byrow = TRUE)) + #split legend in two columns
  theme_bw(base_size=20) + theme(legend.position = "top") #increase font size, put legend on top
dev.off()


##### Plot NA
M <- 10^6
plotNA<-data.frame(cutoff,100*nrNA1/M, 100*nrNA28/M, 100*nrNA55/M, 100*nrNA25/M, 100*nrNA53/M, 100*nrNA79/M) #combine data into one dataframe
colnames(plotNA) <- c("cutoff","k = 2, npilot = 10, ES = 0.0099","k = 2, npilot = 10, ES = 0.0588","k = 2, npilot = 10, ES = 0.1379","k = 4, npilot = 50, ES = 0.0099", "k = 4, npilot = 50, ES = 0.0588", "k = 4, npilot = 50, ES = 0.1379") #2 groups, npilot = 10, small and large ES (and eta; 80% power)
plotNA <- melt(plotNA, id="cutoff") #convert to long format with reshape package to easily plot multiple groups
colnames(plotNA) <- c("cutoff", "condition","value") #change names in final dataframe to plot

tiff(file=paste("FigS1.tiff",sep=""),width=3800,height=2000, units = "px", res = 300)
plot.new()
ggplot(data=plotNA, aes(x = cutoff, y=value, group=condition, colour=condition)) +
  scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(size=1.5) + #geom_point(size=5) + #plot lines and dots
  ylab("% excluded values")  + xlab("cut-off") + #set labels
  #scale_x_continuous(breaks=pretty(plotNA$cutoff, n = 10)) + #use scales package to set breaks
  scale_x_log10(breaks=c(200,500,1000,2000,5000,10000,20000,50000,100000)) + 
  guides(col = guide_legend(nrow = 2, byrow = TRUE)) + #split legend in two columns
  theme_bw(base_size=20) + theme(legend.position = "top") #increase font size, put legend on top
dev.off()


