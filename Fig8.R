library(pwr)
library(reshape2)
library(ggplot2)

min <- 25
max <- 250

nstar <- seq(max,min)  # note: the values *must* be listed in decreasing order

unbiased_power_small<-pwr.anova.test(f=0.1,k=2,n=seq(min,max),sig.level=0.05)$power
unbiased_power_medium<-pwr.anova.test(f=0.25,k=2,n=seq(min,max),sig.level=0.05)$power
unbiased_power_large<-pwr.anova.test(f=0.4,k=2,n=seq(min,max),sig.level=0.05)$power

unbiased_power <- matrix(NA, nrow=3,ncol=length(nstar))
unbiased_power[1,]<-unbiased_power_small
unbiased_power[2,]<-unbiased_power_medium
unbiased_power[3,]<-unbiased_power_large

unbiased_power<-t(unbiased_power)
colnames(unbiased_power) <- c("ES = 0.0099","ES = 0.0588","ES = 0.1379") 

unbiased_power <- melt(unbiased_power, id="ES") #convert to long format with reshape package to easily plot multiple groups
colnames(unbiased_power) <- c("N","ES","Power") 
cbbPalette<-c("#5e3c99", "#e66101", "#fdb863", "#000000")

tiff(file="Fig_8.tiff",width=2700,height=1700, units = "px", res = 300)
ggplot(data=unbiased_power, aes(x = N, y=Power, group=ES)) +
  #scale_color_manual(values=cbbPalette) + #use custum color palette
  geom_line(aes(linetype=ES),size=1.5) + #geom_point(size=5) + #plot lines and dots
  scale_x_continuous(breaks = c(25,75,125,175,225), labels= c(50,100,150,200,250)) + 
  scale_linetype_manual(values=c(1,2,3)) + 
  ylab("Power")  + xlab("Sample size per condition") + #set labels
  theme_bw(base_size=20) + expand_limits(y=c(0,1)) + theme(legend.title=element_blank(), legend.position = "top", legend.key.width=unit(2,"cm")) #increase font size, put legend on top
dev.off()
