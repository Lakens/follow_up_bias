library(ggplot2)
library(pwr)
# Self-built functions

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
    epsilon2 <- max((SSb - (K-1)*MSw)/SSt,0)     # the 'max' makes sure it's non-negative
    omega2 <- max((SSb - (K-1)*MSw)/(SSt+MSw),0)
    ESMAT[r,] <- c(eta2, epsilon2, omega2)
    
}
DichoMAT <- ESMAT < ES
apply(DichoMAT,2,mean)
library("reshape2")
ESs <- melt(ESMAT)
colnames(ESs) <- c("R","ES","value")
ESs <- ESs[ESs$R < 1000,] 
#save(ESMAT,file=paste('ESMAT.RDATA'))
#load(file=paste('ESMAT.RDATA'))
plot.new()
ggplot(data=ESs, aes(x=as.factor(ES), y=value)) + 
  geom_violin(trim=F, fill="red") +
  scale_color_manual(values=cbbPalette) +
  ylab("effect size") +
  xlab("") + 
  theme_bw(base_size=20) + theme(legend.position = "top") 
dev.off()
