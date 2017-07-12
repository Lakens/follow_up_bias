#Create Dataset
muvec<-c(-0.2499469,0.2499469)
sd<-1
cK<-2
n<-50
nsims<-100000

# create progress bar
pb <- winProgressBar(title = "progress bar", min = 0, max = nsims, width = 300)

ESMAT<-numeric(nsims)
for(i in 1:nsims){
  setWinProgressBar(pb, i, title=paste(round(i/(nsims)*100, 2), "% done"))
  study.data  <- rnorm(n = cK*n, mean = rep(muvec,each=n), sd = rep(1,cK*n))
  study.anova <- as.matrix(anova(aov(study.data ~ rep(as.factor(1:cK),each=n))))
  SSb <- study.anova[1,2]
  SSt <- study.anova[1,2] + study.anova[2,2]
  MSw <- study.anova[2,3]
  eta2 <- SSb/SSt
  ESMAT[i] <- eta2 
}
#close progress bar
close(pb)

library(pwr)

ESlow<-ESMAT[ESMAT<mean(ESMAT)]
length(ESlow)/nsims #percentage of effect size estimate below the mean of all effect sizes.
ESlow<-ESMAT[ESMAT<0.0588]
length(ESlow)/nsims #percentage of effect size estimate below the mean of all effect sizes.


etas<-seq(0.025,0.3,0.025)
f2 <- etas/(1-etas)
f<-sqrt(f2)

#f2<-seq(0.1,0.5,0.1)
n_power<-numeric(length(f))
for(i in 1:length(f)){
  n_power[i]<-round(pwr.anova.test(k = 2, f = f[i], sig.level = 0.05, power = 0.8)$n, digits=0)
}



tiff(file=paste("Fig1.tiff",sep=""),width=2000,height=1000, units = "px", res = 300)
par(mar=c(8, 1, 1, 1) + 0.1) #Define Margins. 
hist(ESMAT, axes=F, xlab="", ylab="", main="",xlim=c(0,0.3), breaks=50, col="grey") #histogram
axis(1, xlim=c(0,max(ESMAT)),col="black",lwd=2) #draw first axis
mtext(1,text=expression(eta^2),line=2.4) #drawsecond axis
par(new=T)
axis(1,lwd=2,line=3.7,at=seq(0.025,0.3,0.025),labels=n_power)
mtext(1,text="n (per group) for 80% power",line=5.8)
dev.off()