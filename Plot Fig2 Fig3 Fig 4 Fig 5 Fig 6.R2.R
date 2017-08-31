###########################################################################
# Taking into account bias and variability of effect sizes in power studies
# Lakens & Albers
###########################################################################

rm(list=ls()) #remove all items (to reduce memory load)
gc() #garbage cleaning

# The script in this file analyses the simulation results
setwd("C:\\Users\\DLakens\\surfdrive\\Word\\ESpaper Casper\\Power Casper")

options(scipen=20) #disable scientific notation

i=13 #set to conditions you want to plot. Use 13 for plot 2, 3, and 4, 1 for plot 5, and 7 for plot 6 
#See power_script.R for the script performing the simulations

K      <- c(2, 3, 4)        # number of groups in a one-way anova
ES     <- c(.0099,.0588,.1379) # population ES
npilot <- c(10, 25, 50)     # sample size per group in the pilot
alpha  <- 0.05              # nominal level of significance
beta   <- c(0.2, 0.1)       # 1 - power
ESmeasure <- c("eta","epsilon","omega")
R      <- 1000000              # number of replications per condition. Set this to 10^6 at the end
conditions <- expand.grid(ESmeasure,K, npilot, ES, alpha, 1-beta)
colnames(conditions) <- c("ESmeasure", "K", "npilot", "ES", "alpha", "power")

npilot <- c(10, 25, 50)     # sample size per group in the pilot

if(i == 1){
  powerselection<-0.8
  ESselection<-0.0099
  Kselection<-2
}
if(i == 2){
  powerselection<-0.8
  ESselection<-0.0099
  Kselection<-3
}
if(i == 3){
  powerselection<-0.8
  ESselection<-0.0099
  Kselection<-4
}
if(i == 4){
  powerselection<-0.9
  ESselection<-0.0099
  Kselection<-2
}
if(i == 5){
  powerselection<-0.9
  ESselection<-0.0099
  Kselection<-3
}
if(i == 6){
  powerselection<-0.9
  ESselection<-0.0099
  Kselection<-4
}
if(i == 7){
  powerselection<-0.8
  ESselection<-0.0588
  Kselection<-2
}
if(i == 8){
  powerselection<-0.8
  ESselection<-0.0588
  Kselection<-3
}
if(i == 9){
  powerselection<-0.8
  ESselection<-0.0588
  Kselection<-4
}
if(i == 10){
  powerselection<-0.9
  ESselection<-0.0588
  Kselection<-2
}
if(i == 11){
  powerselection<-0.9
  ESselection<-0.0588
  Kselection<-3
}
if(i == 12){
  powerselection<-0.9
  ESselection<-0.0588
  Kselection<-4
}
if(i == 13){
  powerselection<-0.8
  ESselection<-0.1379
  Kselection<-2
}
if(i == 14){
  powerselection<-0.8
  ESselection<-0.1379
  Kselection<-3
}
if(i == 15){
  powerselection<-0.8
  ESselection<-0.1379
  Kselection<-4
}
if(i == 16){
  powerselection<-0.9
  ESselection<-0.1379
  Kselection<-2
}
if(i == 17){
  powerselection<-0.9
  ESselection<-0.1379
  Kselection<-3
}
if(i == 18){
  powerselection<-0.9
  ESselection<-0.1379
  Kselection<-4
}


# Selecting subset of data (cannot load everything in memory at the same time)
nrconditions=162
conditionselection<-c(1:162) #Select conditions
conditions$conditionselect <- conditionselection #Add numbers for 162 conditions to select them below

conditionselection<-subset(conditions, power==powerselection & ES==ESselection & K==Kselection) #Keep only those conditions from 162 we want to plot
conditionselection<-conditionselection$conditionselect
nrconditions=length(conditionselection) #get number of condition we are plotting

conditions<-conditions[conditionselection,] #remove conditions we are not interested in 

# 4. Load the conditions you want to study ----
nMAT <- array(NA, dim=c(nrconditions,R))
ESMAT <- nMAT
pvalueMAT <- nMAT 
for(i in 1:length(conditionselection)){
  cond<-conditionselection[i]
  nMAT[i,] <- readRDS(paste("n_",cond,".rds",sep=""))
  ESMAT[i,] <- readRDS(paste("ES_",cond,".rds",sep=""))
  pvalueMAT[i,] <- readRDS(paste("p_",cond,".rds",sep=""))
}


# 5. Some summaries over replications ----
# The 'raw' simulation results are retained as well, e.g. in nMAT, which could be used for e.g. a violin plot
nmedian2 <- apply(nMAT,1,median, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nmean2 <- apply(nMAT,1,mean, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
large_n <- rowSums(is.na(nMAT))/R   # percentage of values with n > 100,000
nnMAT <- nMAT
nnMAT[is.na(nMAT)] <- 10^6
rm(nMAT) #clean up to save memory

nmedian <- apply(nnMAT,1,median, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nmean <- apply(nnMAT,1,mean, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
biasES <- apply(ESMAT,1,mean, na.rm=TRUE) - conditions$ES
meanES <- apply(ESMAT,1,mean, na.rm=TRUE)
sdES <- apply(ESMAT,1,sd, na.rm=TRUE)
rmseES <- sqrt(biasES^2 + sdES^2)
# Power is now calculated by dividing by R - but we should divide by R - NA's: obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)
obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/R*100) #kept, but not used, incorrect power
obspower2 <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x))))
results <- cbind(conditions,nmedian,nmedian2,nmean,nmean2,biasES,meanES,sdES,rmseES,large_n, obspower, obspower2)
print(results)

#load plot packages
library(scales)
library(ggplot2)
library(gridExtra)
library(pwr)

#create color scheme for plots
cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #old color scheme
cbbPalette<-c("#5e3c99", "#e66101", "#fdb863", "#000000") # new print and colorblind safe color scheme

#Create dataframe for plots, specify all factors
names<-paste("Groups =", results$K, "Pilot n =", results$npilot, "ES =", results$ES, "Power =", results$power, results$ESmeasure) #Create name labels

nameslist<-rep(names, each = R)
npilotlist<-rep(results$npilot, each = R)
ESconlist<-rep(results$ESmeasure, each = R)
powerlist<-rep(results$power, each = R)
print(names)
#add data: ES, n, p-value
ESlist<-as.vector(t(ESMAT))
nlist<-as.vector(t(nnMAT))
plist<-as.vector(t(pvalueMAT))
plotdata<-data.frame(nameslist, ESlist, nlist, plist, npilotlist, ESconlist, powerlist)
rm(nameslist, ESlist, nlist, plist, npilotlist, ESconlist, powerlist) #clean up to save memory
colnames(plotdata) <- c("condition", "ES", "SS", "pvalue", "NPilot", "EffectSize", "Power")

results <- cbind(names,conditions,nmedian,nmedian2,nmean,nmean2,biasES,meanES,sdES,rmseES,large_n, obspower, obspower2)

plotdata<-na.omit(plotdata) #remove NA


#Plot effect sizes ----
tiff(file=paste("Fig2_",toString(conditionselection),".tiff",sep=""),width=2500,height=1700, units = "px", res = 300)
plot.new()
eval(bquote(ggplot(data=plotdata, aes(factor(condition),ES, fill = factor(plotdata$EffectSize))))) + 
  geom_boxplot(outlier.colour=alpha("black", 0.2)) + scale_fill_manual(values=cbbPalette, name = "Effect Size", labels=c(parse(text="   epsilon^2"), parse(text="   eta^2"), parse(text="   omega^2")), breaks=c("epsilon","eta","omega")) +
  facet_grid(.~NPilot, scales="free", space="free") +
  geom_violin(alpha=0.5) +
  theme(text=element_text(size=18), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        axis.title.x=element_blank(), strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white'), legend.text.align = 0, legend.text = element_text(size = rel(1.0))) +
  geom_hline(yintercept=unique(results$ES), colour="gray20", linetype="dashed",size=1) +
  stat_summary(fun.y=mean, colour="#000000", geom="point", size=3) +
  ylab("Effect Size Estimate")
dev.off()


#Plot sample sizes ----
#First remove high N above cut-off
cutoff <- 250
plotdata$SS[plotdata$SS>cutoff] <- NA
percent <- round(sum(is.na(plotdata$SS))/length(plotdata$SS),3)*100
labelmissing <- paste(percent,"% over 250",sep="")
plotdata<-subset(plotdata, !is.na(plotdata$SS))

cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #old color scheme

#plot sample sizes ----
tiff(file=paste("Fig3_",toString(conditionselection),".tiff",sep=""),width=2500,height=1700, units = "px", res = 300)
plot.new()
ggplot(data=plotdata, aes(factor(condition),SS, fill = factor(EffectSize))) + 
  #geom_boxplot(outlier.colour=alpha("white", 0.05)) + 
  scale_fill_manual(values=cbbPalette, name = "Effect Size", labels=c(parse(text="epsilon^2"), parse(text="eta^2"), parse(text="omega^2")), breaks=c("epsilon","eta","omega")) +
  facet_grid(.~NPilot, scales="free", space="free") +
  geom_violin(alpha=0.95) +
  theme(text=element_text(size=18), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), 
        axis.title.x=element_blank(), panel.background = element_rect(fill = "white"), strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white'), legend.text.align = 0, legend.text = element_text(size = rel(1.0))) +
  geom_hline(yintercept=ceiling(pwr.anova.test(f=sqrt(ESselection/(1-ESselection)),k=Kselection,sig.level=0.05, power=powerselection)$n), colour="gray20", linetype="dashed",size=1) +
  stat_summary(fun.y=mean, colour="#000000", geom="point", size=3, shape = 18) +
  ylab("Sample Size")
dev.off()

cbbPalette<-c("#5e3c99", "#e66101", "#fdb863", "#000000") # new print and colorblind safe color scheme


#Plot observed power----
tiff(file=paste("Fig4_",toString(conditionselection),".tiff",sep=""),width=2500,height=1700, units = "px", res = 300)
plot.new()
ggplot(results, aes(x=names, y=obspower2, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size", labels=c(parse(text="epsilon^2"), parse(text="eta^2"), parse(text="omega^2")), breaks=c("epsilon","eta","omega")) +
  facet_grid(.~npilot, scales="free", space="free") +
  theme(text=element_text(size=18), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        axis.title.x=element_blank(), strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white'), legend.text.align = 0, legend.text = element_text(size = rel(1.0))) +
  geom_hline(yintercept=(unique(results$power)), colour="gray20", linetype="dashed", size = 1.4) +
  ylab("Observed Power")  
dev.off()
####