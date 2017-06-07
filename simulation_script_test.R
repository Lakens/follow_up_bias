###########################################################################
# Taking into account bias and variability of effect sizes in power studies
# Lakens & Albers
###########################################################################
# Parts of the code have been recycled from 
# Okada, K. (2013), Behaviormetrika, 40(2), 129-147.


# TABLE OF CONTENTS
# 0. Prepare environment
# 1. Description of the simulation study
# 2. Setting up parameters of the study
# 3. Running the simulations

# 0. Prepare environment ----
# Load required packages:
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


# 1. Description of the simulation study ----
# (Here we can add a few lines of explanation / the abstract of the paper / in comments)

# 2. Setting up parameters of the study ----

K      <- c(2, 3, 4)        # number of groups in a one-way anova
npilot <- c(10, 25, 50)     # sample size per group in the pilot
ES     <- c(.0099,.0588,.1379) # population ES
alpha  <- 0.05              # nominal level of significance
beta   <- c(0.2, 0.1)       # 1 - power
ESmeasure <- c("eta","epsilon","omega")
R      <- 1              # number of replications per condition. Set this to 10^6 at the end
conditions <- expand.grid(ESmeasure,K, npilot, ES, alpha, 1-beta)
colnames(conditions) <- c("ESmeasure", "K", "npilot", "ES", "alpha", "power")
nrconditions <- dim(conditions)[1]
#for(i in 1:nrconditions){
#  conditions[i,"muvec"] <- selectmuvec(conditions[i,"K"], conditions[i,"ES"])
#}


# 3. Run the simulations ----
simulateconditions <- 37:37 # change 1:nrconditions in the range you want to simulate on your machine; e.g. 11:25
for(cond in simulateconditions){ 
#  set.seed(cond)     # set the seed, such that the study's results can be reproduced
  nMAT <- rep(NA, R)
  ESMAT <- rep(NA, R)
  pvalueMAT <- rep(NA, R)
  
  muvec <- selectmuvec(conditions[cond,"K"], conditions[cond,"ES"])
  cK <- conditions[cond,"K"] # abbrevated notation, because this one's used often
  cnp <- conditions[cond,"npilot"] # same
  cES <- conditions[cond,"ES"] # same
  for(r in 1:R){
    # step 1: perform the pilot and estimate ES
    pilot.data  <- rnorm(n = cK*cnp, mean = rep(muvec,each=cnp), sd = rep(1,cK*cnp))
    pilot.anova <- as.matrix(anova(aov(pilot.data ~ rep(as.factor(1:cK),each=cnp))))
    SSb <- pilot.anova[1,2]
    SSt <- pilot.anova[1,2] + pilot.anova[2,2]
    MSw <- pilot.anova[2,3]
    eta2p <- SSb/SSt
    epsilon2 <- max((SSb - (cK-1)*MSw)/SSt,0)     # the 'max' makes sure it's non-negative
    omega2 <- max((SSb - (cK-1)*MSw)/(SSt+MSw),0)
    
    # step 2: compute the required sample size
    if(conditions[cond, "ESmeasure"] == "eta"){
      tmp <-  try(pwr.anova.test(f=sqrt(eta2/(1-eta2)),k=cK,sig.level=0.05, power=conditions[cond,"power"])$n, silent = TRUE)
      if ('try-error' %in% class(tmp)){
        n <- NA
      }else{
        n <- ceiling(tmp) # round it to an integer. Round it up to ensure that you achieve the power
      }
    }
    
    if(conditions[cond, "ESmeasure"] == "epsilon"){
      tmp <-  try(pwr.anova.test(f=sqrt(epsilon2/(1-epsilon2)),k=cK,sig.level=0.05, power=conditions[cond,"power"])$n, silent = TRUE)
      if ('try-error' %in% class(tmp)){
        n <- NA
      }else{
        n <- ceiling(tmp) # round it to an integer. Round it up to ensure that you achieve the power
      }
    }
  
    if(conditions[cond, "ESmeasure"] == "omega"){
      tmp <-  try(pwr.anova.test(f=sqrt(omega2/(1-omega2)),k=cK,sig.level=0.05, power=conditions[cond,"power"])$n, silent = TRUE)
      if ('try-error' %in% class(tmp)){
        n <- NA
      }else{
        n <- ceiling(tmp) # round it to an integer. Round it up to ensure that you achieve the power
      }
    }

    # Step 3: perform the analysis with req. sample size, then compute an ES-estimate
    if(is.na(n)){
      eta2 <- epsilon2 <- omega2 <- NA
    }else{
      study.data  <- rnorm(n = cK*n, mean = rep(muvec,each=n), sd = rep(1,cK*n))
      study.anova <- as.matrix(anova(aov(study.data ~ rep(as.factor(1:cK),each=n))))
      SSb <- study.anova[1,2]
      SSt <- study.anova[1,2] + study.anova[2,2]
      MSw <- study.anova[2,3]
      eta2 <- SSb/SSt
      epsilon2 <- (SSb - (cK-1)*MSw)/SSt
      omega2 <- (SSb - (cK-1)*MSw)/(SSt+MSw)
      pvalueMAT[r] <- study.anova[1,5]
    }
    # Step 4: save what we need to save
    nMAT[r] <- n
    if(conditions[cond, "ESmeasure"] == "eta"){ 
      ESMAT[r] <- eta2 
    }
    if(conditions[cond, "ESmeasure"] == "epsilon"){ 
      ESMAT[r] <- epsilon2 
    }
    if(conditions[cond, "ESmeasure"] == "omega"){ 
      ESMAT[r] <- omega2 
    }
  }
  # Save results and two images
  #saveRDS(nMAT,paste("n_",cond,".rds",sep=""))   # Save each simulation-condition to a separate file, for easy parallel simulation
  #saveRDS(ESMAT,paste("ES_",cond,".rds",sep=""))
  #saveRDS(pvalueMAT, paste("p_",cond,".rds",sep=""))
}

# The script above performs all simulations
######################################################################################
eta2p
eta2