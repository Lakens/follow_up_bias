###########################################################################
# Taking into account bias and variability of effect sizes in power studies
# Lakens & Albers
###########################################################################

rm(list=ls())
gc()

options(scipen=20) #disable scientific notation

#See power_script.R for the script performing the simulations

K      <- c(2, 3, 4)        # number of groups in a one-way anova
npilot <- c(10, 25, 50)     # sample size per group in the pilot
ES     <- c(.0099,.0588,.1379) # population ES
alpha  <- 0.05              # nominal level of significance
beta   <- c(0.2, 0.1)       # 1 - power
ESmeasure <- c("eta","epsilon","omega")
R      <- 1000000              # number of replications per condition. Set this to 10^6 at the end
conditions <- expand.grid(ESmeasure,K, npilot, ES, alpha, 1-beta)
colnames(conditions) <- c("ESmeasure", "K", "npilot", "ES", "alpha", "power")

# The script in this file, analyses the simulation results
setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")

# Selecting subset of data (cannot load everything in memory at the same time)
R=1000000
nrconditions=162
conditionselection<-c(1:162) #Select conditions

# 4. Load the conditions you want to study ----
nMAT <- array(NA, dim=c(nrconditions,R))
ESMAT <- array(NA, dim=c(nrconditions,R))
pvalueMAT <- array(NA, dim=c(nrconditions,R)) 
for(cond in conditionselection){
  nMAT[cond,] <- readRDS(paste("n_",cond,".rds",sep=""))
  ESMAT[cond,] <- readRDS(paste("ES_",cond,".rds",sep=""))
  pvalueMAT[cond,] <- readRDS(paste("p_",cond,".rds",sep=""))
}

# Count missing values in each condition

CountNAn<-apply(nMAT,1,function(x) sum(is.na(x)))
CountNAES<-apply(ESMAT,1,function(x) sum(is.na(x)))
CountNAp<-apply(pvalueMAT,1,function(x) sum(is.na(x)))
CountNA<-cbind(conditions,CountNAn,CountNAES,CountNAp)
write.table(CountNA, "countNA.txt")

#Count maximum sample size
max(nMAT, na.rm=TRUE)
mean(nMAT, na.rm=TRUE)
median(nMAT, na.rm=TRUE)
table(nMAT)
hist(nMAT, breaks=10000)
hist(nMAT, breaks=1000, ylim=c(0, 10000))
hist(nMAT, breaks=10000, xlim=c(0, 1000))



#plot observed effect sizes
ESdataplotES<-c(ESMAT[1,],ESMAT[2,],ESMAT[3,])
ESdataplotnames<-as.factor(rep(c("eta", "epsilon", "omega"), each = 1000000))
ESdataplot<-as.data.frame(cbind(ESdataplotnames,ESdataplotES))
head(ESdataplot)

library(ggplot2)
ggplot(ESdataplot, aes(x = ESdataplotES, na.rm = TRUE)) + 
  geom_histogram(binwidth = 0.01) + 
  theme(text=element_text(size=16), panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white'))   + coord_cartesian(xlim=c(-0.2,0.6)) +
  facet_grid(.~ESdataplotnames, scales="free", space="free")


#plot subset of observed effect sizes cut-off by SESOI----
ESmeans<-apply(ESMAT,1,mean, na.rm=TRUE)
SESOI<-0
ESMAT <- ifelse(ESMAT<SESOI,NA,ESMAT)
ESmeansCUT<-apply(ESMAT,1,mean, na.rm=TRUE)
ESmeansTOTAL<-as.data.frame(cbind(ESmeans,ESmeansCUT))
names<-paste("Groups =", conditions$K, "Pilot n =", conditions$npilot, "ES =", conditions$ES, "Power =", conditions$power, conditions$ESmeasure) #Create name labels
ESmeansTOTAL$names<-names
ESmeansTOTAL$Dif<-ESmeansTOTAL$ESmeansCUT-ESmeansTOTAL$ESmeans



library(scales)
library(ggplot2)
library(gridExtra)
#Plot ALL observed differences in ES----
png(file=paste("obsDifESALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(ESmeansTOTAL, aes(x=names, y=Dif)) + 
  geom_bar(stat="identity") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  ylab("Observed Power")  + xlab("Condition")  + coord_cartesian(ylim=c(0,0.1)) + #set labels and max of y-axis
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))
p3
dev.off()
####
####

ESdataplotES<-c(ESMAT[1,],ESMAT[2,],ESMAT[3,])
ESdataplotnames<-as.factor(rep(c("eta", "epsilon", "omega"), each = 1000000))
ESdataplot<-as.data.frame(cbind(ESdataplotnames,ESdataplotES))

library(ggplot2)
ggplot(ESdataplot, aes(x = ESdataplotES, na.rm = TRUE)) + 
  geom_histogram(binwidth = 0.01) + 
  theme(text=element_text(size=16), panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white'))   + coord_cartesian(xlim=c(-0.2,0.6)) +
  facet_grid(ESdataplotnames~.)



##Data unavailable, this could be run on original ES
#ESdataplotDist<-c(subset(ESMAT[1,],ESMAT[1,]>0))
#mean(ESdataplotDist)
#mean(ESMAT[1,],na.rm=TRUE)

# 5. Some summaries over replications ----
# The 'raw' simulation results are retained as well, e.g. in nMAT, which could be used for e.g. a violin plot
large_n <- rowSums(is.na(nMAT))/R   # percentage of values with n > 100,000
nnMAT <- nMAT
nnMAT[is.na(nMAT)] <- 10^6
rm(nMAT) #clean up to save memory
nmedian <- apply(nnMAT,1,median, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nQ1 <- apply(nnMAT,1,quantile, probs = .25, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nQ3 <- apply(nnMAT,1,quantile, probs = .75, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nmean <- apply(nnMAT,1,mean, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
biasES <- apply(ESMAT,1,mean, na.rm=TRUE) - conditions$ES
sdES <- apply(ESMAT,1,sd, na.rm=TRUE)
rmseES <- sqrt(biasES^2 + sdES^2)
# Power is now calculated by dividing by R - but we should divide by R - NA's: obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)
obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/R*100)
obspower2 <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)
results <- cbind(conditions,nmedian,nmean,biasES,sdES,rmseES,large_n, obspower, obspower2)
names<-paste("Groups =", results$K, "Pilot n =", results$npilot, "ES =", results$ES, "Power =", results$power, results$ESmeasure) #Create name labels
results <- cbind(names,conditions,nmedian,nmean,biasES,sdES,rmseES,large_n, obspower, obspower2)

print(results)
write.table(results, "results.txt")


# 6. Describing the results ----

#load plot packages
library(scales)
library(ggplot2)
library(gridExtra)

#create color scheme for plots
cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



#Plot ALL observed power----
png(file=paste("obspowerALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(results, aes(x=names, y=obspower2, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  geom_hline(yintercept=(100*unique(results$power)), colour="gray20", linetype="dashed") +
  ylab("Observed Power")  + xlab("Condition")  + coord_cartesian(ylim=c(0,100)) + #set labels and max of y-axis
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1)) + scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####

#Plot ALL observed means N----
png(file=paste("obsnmeanALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(results, aes(x=names, y=nmean, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  ylab("Mean N")  + xlab("Condition")  + coord_cartesian(ylim=c(0,250)) + #set labels and max of y-axis
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1)) + scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####


#Plot ALL observed median N----
png(file=paste("obsnmedianALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(results, aes(x=names, y=nmedian, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  ylab("Median N")  + xlab("Condition")  + coord_cartesian(ylim=c(0,250)) + #set labels and max of y-axis
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1)) + scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####

#Plot ALL observed bias ES----
png(file=paste("obsbiasESALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(results, aes(x=names, y=biasES, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  ylab("Bias ES")  + xlab("Condition")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1)) + scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####


#Plot ALL observed RMSE----
png(file=paste("obsRMSE_ALL.png",sep=""),width=9000,height=3000, res = 300)
plot.new()
p3 <- ggplot(results, aes(x=names, y=rmseES, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white")) +
  ylab("RMSE")  + xlab("Condition")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1)) + scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####

#check correlations between output of the simulation
cor_results<-cbind(nmedian,nmean,biasES,sdES,rmseES,large_n, obspower, obspower2)
cor(cor_results)

