###########################################################################
# Taking into account bias and variability of effect sizes in power studies
# Lakens & Albers
###########################################################################

rm(list=ls())
gc()
i=1 #Poor coders loop - can't get the real loop to work, so run it 18 times changing this from 1:18
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
R=1000000
nrconditions=162
conditionselection<-c(1:162) #Select conditions
conditions$conditionselect <- conditionselection

conditionselection<-subset(conditions, power==powerselection & ES==ESselection & K==Kselection)
conditionselection<-conditionselection$conditionselect

# 4. Load the conditions you want to study ----
nMAT <- array(NA, dim=c(nrconditions,R))
ESMAT <- nMAT
pvalueMAT <- nMAT 
for(cond in conditionselection){
  nMAT[cond,] <- readRDS(paste("n_",cond,".rds",sep=""))
  ESMAT[cond,] <- readRDS(paste("ES_",cond,".rds",sep=""))
  pvalueMAT[cond,] <- readRDS(paste("p_",cond,".rds",sep=""))
}

#save only selected conditions to reduce memory use ----
nMAT<-nMAT[conditionselection,]
ESMAT<-ESMAT[conditionselection,]
pvalueMAT<-pvalueMAT[conditionselection,]
conditions<-conditions[conditionselection,]


# 5. Some summaries over replications ----
# The 'raw' simulation results are retained as well, e.g. in nMAT, which could be used for e.g. a violin plot
large_n <- rowSums(is.na(nMAT))/R   # percentage of values with n > 100,000
nnMAT <- nMAT
nnMAT[is.na(nMAT)] <- 10^6
#rm(nMAT) #clean up to save memory
nmedian <- apply(nnMAT,1,median, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nQ1 <- apply(nnMAT,1,quantile, probs = .25, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
nQ3 <- apply(nnMAT,1,quantile, probs = .75, na.rm=TRUE)  # I prefer median over mean, since distribution is so skewed
biasES <- apply(ESMAT,1,mean, na.rm=TRUE) - conditions$ES
sdES <- apply(ESMAT,1,sd, na.rm=TRUE)
rmseES <- sqrt(biasES^2 + sdES^2)
# Power is now calculated by dividing by R - but we should divide by R - NA's: obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)
obspower <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/R*100)
obspower2 <- apply(pvalueMAT,1, function(x) sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)
results <- cbind(conditions,nmedian,biasES,sdES,rmseES,large_n, obspower, obspower2)
print(results)
write.table(results, "results.txt")
# 6. Describing the results ----
options(scipen=20)
apply(pvalueMAT, 1, function(x) hist(x, main = paste(sum(x < .05, na.rm=TRUE)/(R-sum(is.na(x)))*100)))

#load plot packages
library(scales)
library(ggplot2)
library(gridExtra)

#create color scheme for plots
cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
results <- cbind(names,conditions,nmedian,biasES,sdES,rmseES,large_n, obspower, obspower2)

#Plot effect sizes ----
png(file=paste("ESestimate_",toString(conditionselection),".png",sep=""),width=1800,height=1250, res = 150)
plot.new()
p1 <- eval(bquote(ggplot(data=plotdata, aes(factor(condition),ES, fill = factor(plotdata$EffectSize))))) + 
  geom_boxplot(outlier.colour=alpha("black", 0.2)) + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  facet_grid(.~NPilot, scales="free", space="free") +
  geom_violin(alpha=0.5) +
  theme(text=element_text(size=16), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white')) +
  geom_hline(yintercept=unique(results$ES), colour="gray20", linetype="dashed",size=1) +
  stat_summary(fun.y=mean, colour="#000000", geom="point", 
               size=3,show_guide = FALSE) +
  ylab("Effect Size")  + xlab(paste("Groups = ",Kselection,", ES = ",ESselection,", Power = ",powerselection,sep=""))
p1
dev.off()

#Plot sample sizes ----
#First remove high N above cut-off
cutoff <- 250
plotdata$SS[plotdata$SS>cutoff] <- NA
percent <- round(sum(is.na(plotdata$SS))/length(plotdata$SS),3)*100
labelmissing <- paste(percent,"% over 250",sep="")
plotdata<-subset(plotdata, !is.na(plotdata$SS))

#plot sample sizes ----
png(file=paste("Nestimate_",toString(conditionselection),".png",sep=""),width=1800,height=1250, res = 150)
plot.new()
p2 <- ggplot(data=plotdata, aes(factor(condition),SS, fill = factor(EffectSize))) + 
  geom_boxplot(outlier.colour=alpha("black", 0.2)) + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  facet_grid(.~NPilot, scales="free", space="free") +
  geom_violin(alpha=0.5) +
  theme(text=element_text(size=16), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), 
        panel.background = element_rect(fill = "white"), strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white')) +
  #  geom_hline(yintercept=unique(results$ES), colour="gray20", linetype="dashed",size=1) +
  stat_summary(fun.y=mean, colour="#000000", geom="point", 
               size=3,show_guide = FALSE) +
  ylab("Sample Size")  + xlab(paste("Groups = ",Kselection,", ES = ",ESselection,", Power = ",powerselection, ", ",labelmissing, sep=""))
p2
dev.off()


#Plot observed power----
png(file=paste("obspower_",toString(conditionselection),".png",sep=""),width=1800,height=1250, res = 150)
plot.new()
p3 <- ggplot(results, aes(x=names, y=obspower2, group=1, fill = factor(ESmeasure))) + 
  geom_bar(stat="identity") + scale_fill_manual(values=cbbPalette, name = "Effect Size") +
  facet_grid(.~npilot, scales="free", space="free") +
  theme(text=element_text(size=16), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16), strip.background = element_rect(fill = 'white')) +
  geom_hline(yintercept=(100*unique(results$power)), colour="gray20", linetype="dashed") +
  ylab("Observed Power")  + xlab(paste("Groups = ",Kselection,", ES = ",ESselection,", Power = ",powerselection,sep=""))  + coord_cartesian(ylim=c(0,100)) + #set labels and max of y-axis
  scale_y_continuous(breaks=pretty_breaks(n=10))
p3
dev.off()
####
