################################################################################
#          Simulation and Model fitting for the PCFG multievent model          #
#           Initial State Probability a Function of Minimum Tenure             #
# Based on code from O. Gimenez https://github.com/oliviergimenez/banana-book/ #
################################################################################


#====================================
# Parameter and Population Settings
#====================================

# Population and Time -----------------
N      <- 300
time   <- 10
states <- c(1, 2, 3, 4) # 1-calf, 2-PCFG, 3-ENP, 4-dead/transient
obs    <- c(1, 2, 3)    # 1-ascertained calf, 2-ascertained PCFG, 3-unascertained


# Parameters -------------------------------
# Initial state probabilities
# Prob of being PCFG function of avg residence time
# Prob of being ENP is 1 - pPCFG
# Prob of calf function of expected fecundity
init.C <- 0.24     # 0.24/1.24 = 0.20

# Initial observation probability, depends on state
#   Calves assumed to be observed correctly
#   All adults are unascertained at first observation

# Transition probabilities
ps <- matrix(c(0.0, 0.7, 0.0, 0.3,
               0.0, 0.9, 0.0, 0.1,
               0.0, 0.2, 0.0, 0.8,
               0.0, 0.0, 0.0, 1.0),
             nrow = length(states), byrow = TRUE)

# Capture probabilities
po <- c(0.7, 0.8, 0.65, 0.0)


#===================================
# Simulation
#===================================
# packages
library(LaplacesDemon)

# Set-up data sets
f  <- sample(1:9, N, replace=TRUE)               # first observation in data set
MT <- sample(1:182, N, replace=T)               # Minimum tenure
Z  <- data.frame(matrix(NA,nrow=N,ncol=time))  # true underlying state
Y  <- data.frame(matrix(0,nrow=N,ncol=time))   # observation


for (i in 1:N){
  # Define latent state at first capture
  init.P <- invlogit(-5 + 0.05*MT[i])
  init.E <- 1 - init.P
  initial_prob <- c(init.C, init.P, init.E)/ (init.C + init.P + init.E)
  Z[i,f[i]] <- rcat(1, initial_prob)
  Y[i,f[i]] <- ifelse(Z[i,f[i]]==1, 1, 3)
  for (t in (f[i]+1):time){
    # State process: draw S(t) given S(t-1)
    Z[i,t] <- rcat(1, ps[Z[i,t-1], 1:4])
    # Observation process: draw O(t) given S(t)
    Obs   <- rbinom(1, 1, po[Z[i,t]])
    if(Z[i,t]==0){
      Zt <- 0
    }else{
      Zt <- ifelse(Z[i,f[i]]==1, 2, 3)
    }
    Y[i,t] <- Zt*Obs
  } #t
} #i



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Model Fitting
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#===================================
# Set up workspace
#===================================

# Packages
library(nimble)
library(MCMCvis)
library(tidyverse)

# Don't print confusing and misleading messages
nimbleOptions(verbose=FALSE)

# Data
y <- Y

# -------------------------------------------------
# PARAMETERS
# pC  detection prob. of calves
# pP  detection prob. of PCFG
# pE  detection prob. of ENP
# psiCP transition prob. from calf to PCFG
# psiEP transition prob. from ENP to PCFG
# psiPD transition prob. from PCFG to transient
# piC prob. of being in initial state Calf
# beta_0 intercept of prob of being PCFG
# beta_1 slope of prob of being PCFG
#-------------------------------------------------
# OBSERVATIONS (y)
# 1 = non-detected
# 2 = seen and ascertained as calf
# 3 = seen and ascertained as PCFG
# 4 = seen and ascertained as ENP (Not possible)
# 5 = not ascertained 
#-------------------------------------------------
#STATES (z)
# 1 = alive calf
# 2 = alive PCFG
# 3 = alive ENP
# 4 = dead/transient
# -------------------------------------------------

unascert.R <- nimbleCode({
  
  # Priors
  psiCP ~ dunif(0, 1)   # Prob transition of calf to PCFG (e.g., survival)
  psiPP ~ dunif(0, 1)   # Prob transition of PCFG to PCFG (e.g., survival)
  psiEP ~ dunif(0, 1)   # Prob transition of ENP to PCFG (e.g., permanent emigration)
  pC ~ dunif(0, 1)      # Prob observing PCFG calf
  pP ~ dunif(0, 1)      # Prob observing PCFG
  pE ~ dunif(0, 1)      # Prob observing ENP
  piC ~ dbeta(1,3)      # Initial state calculation for calf
  beta_0 ~ dnorm(0,1)   # Intercept for prob PCFG model
  beta_1 ~ dnorm(0,1)   # Slope for prob PCFG model
  
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1:4] <- c(0, psiCP, 0, (1 - psiCP))
  gamma[2,1:4] <- c(0, psiPP, 0, (1 - psiPP))
  gamma[3,1:4] <- c(0, psiEP, 0, (1 - psiEP))
  gamma[4,1:4] <- c(0, 0, 0, 1)

  
  # probabilities of y(t) given z(t)
  omega[1,1:5,1] <- c((1 - pC), pC, 0,  0, 0)
  omega[2,1:5,1] <- c((1 - pP), 0,  pP, 0, 0)
  omega[3,1:5,1] <- c((1 - pE), 0,  0 , 0, pE)            
  omega[4,1:5,1] <- c(1, 0, 0, 0, 0)
  omega[1,1:5,2] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,2] <- c((1 - pP), 0,  0, 0, pP)
  omega[3,1:5,2] <- c((1 - pE), 0,  0, 0, pE)            
  omega[4,1:5,2] <- c(1, 0, 0, 0, 0)
  omega[1,1:5,3] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,3] <- c((1 - pP), 0,  0, 0, pP)
  omega[3,1:5,3] <- c((1 - pE), 0,  0, 0, pE)            
  omega[4,1:5,3] <- c(1, 0, 0, 0, 0)
  
  # Likelihood 
  for (i in 1:N){
    # Define initial state probabilities
    piP[i] <- ilogit(beta_0 + beta_1*MT[i])
    piE[i] <- 1 - piP[i]
    delta[1,i] <- piC/(piP[i] + piE[i] + piC)
    delta[2,i] <- piE[i]/(piP[i] + piE[i] + piC)
    delta[3,i] <- piP[i]/(piP[i] + piE[i] + piC)
    delta[4,i] <- 0
    
    # Define latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:4,i])
    y[i,first[i]] <- omega.init[z[i,first[i]]]
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:5,z[i,first[i]]])
    }
  }
})
  


# constants
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# data
yy <- y + 1
yy[yy==4] <- 5
omega.init <- c(1, 5, 5, 0)
my.data <- list(y = yy,
                omega.init = omega.init,
                MT = MT)


# Initial values for Z
zinit <- y

for(i in 1:nrow(y)) {
  for (j in 1:ncol(y)){
    if(j == first[i] & y[i,first[i]]==3) {zinit[i,j] <- sample(c(2,3),1)}
    if(j > first[i]) {zinit[i,j] <- 2}
    if(j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)

# All initial values
initial.values <- function(){list(psiCP = runif(1, 0, 1),
                                  psiPP = runif(1, 0, 1), 
                                  psiEP = runif(1, 0, 1), 
                                  pC = runif(1, 0, 1),
                                  pP = runif(1, 0, 1), 
                                  pE = runif(1, 0, 1),
                                  piC = rbeta(1, 1, 3),
                                  beta_0 = rnorm(1, 0, 1),
                                  beta_1 = rnorm(1, 0, 1),
                                  z = zinit)}  

# MCMC settings
parameters.to.save <- c("psiCP",
                        "psiPP", 
                        "psiEP",
                        "piC",
                        "beta_0",
                        "beta_1",
                        "pC",
                        "pP", 
                        "pE")

n.iter   <- 500000
n.burnin <- 100000
n.chains <- 3
n.thin   <- 10


# Run the models
mcmc.unascertR <- nimbleMCMC(code = unascert.R, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            thin = n.thin,
                            nburnin = n.burnin, 
                            nchains = n.chains)

MCMCsummary(mcmc.unascertR, round = 2)

MCMCplot(mcmc.unascertR)

MCMCtrace(mcmc.unascertR,
          params = c("psiCP",
                       "psiPP", 
                       "psiEP",
                       "piC",
                       "beta_0",
                       "beta_1",
                       "pC",
                       "pP", 
                       "pE"))
