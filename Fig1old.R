library(ggplot2)
library(pwr)
# Self-built functions
setwd("C:/Users/Daniel/BitTorrent Sync/Power Casper")
selectmuvec <- function(K, ES){ # provides the vector of population means for a given population ES and nr of groups
  f2 <- ES/(1-ES)
  if(K == 2){
    a <- sqrt(f2)
    muvec <- c(-a,a)
  }
  if(K == 3){
    a <- sqrt(3*f2/2)
    muvec <- c(-a, 0, a)
  }
  if(K == 4){
    a <- sqrt(f2)
    muvec <- c(-a, -a, a, a)
  } # note: function gives error when K not 2,3,4. But we don't need other K.
  return(muvec)
}



K      <- 3 
#npilot <- c(10, 25, 50)     # sample size per group in the pilot
nmain <- 159/3 # 53 per groep, dat hoort bij verwachte ES medium
ES     <- c(.0588) #,.0588,.1379) # population ES
alpha  <- 0.05              # nominal level of significance
beta   <- c(0.2) #, 0.1)       # 1 - power
ESmeasure <- c("eta","epsilon","omega")
R      <- 1000000   
conditions <- expand.grid(ESmeasure,K, nmain, ES, alpha, 1-beta)
colnames(conditions) <- c("ESmeasure", "K", "n", "ES", "alpha", "power")
nrconditions <- dim(conditions)[1]

muvec <- selectmuvec(K,ES)
ESMAT <- matrix(NA,nrow=R,ncol=3)
set.seed(1)
for(r in 1:R){
    data  <- rnorm(n = K*nmain, mean = rep(muvec,each=nmain), sd = rep(1,K*nmain))
    anova <- as.matrix(anova(aov(data ~ rep(as.factor(1:K),each=nmain))))
    SSb <- anova[1,2]
    SSt <- anova[1,2] + anova[2,2]
    MSw <- anova[2,3]
    eta2 <- SSb/SSt
    epsilon2 <- (SSb - (K-1)*MSw)/SSt     # the 'max' makes sure it's non-negative
    omega2 <- (SSb - (K-1)*MSw)/(SSt+MSw)
    ESMAT[r,] <- c(eta2, epsilon2, omega2)
    
}
DichoMAT <- ESMAT < ES
apply(DichoMAT,2,mean)
library("reshape2")
ESs <- melt(ESMAT)
colnames(ESs) <- c("R","ES","value")
#save(ESMAT,file=paste('ESMAT.RDATA'))
#load(file="ESMAT.RDATA")
ESs$ES <- as.factor(ESs$ES)
#create color scheme for plots
cbbPalette<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(ESs)

tiff(file="Fig1.tiff",width=2500,height=1500, units = "px", res = 300)
plot.new()
ggplot(data=ESs, aes(x=ES, y=value, fill = ES)) + 
  geom_violin() + scale_fill_manual(values=cbbPalette, guide=FALSE) +
  ylab("effect size") +
  xlab("") + scale_x_discrete(limits=c("2","1","3"),labels=c("1" = parse(text="epsilon^2"), "2" = parse(text="eta^2"), "3" = parse(text="omega^2"))) +
  theme_bw(base_size=20)
dev.off()

