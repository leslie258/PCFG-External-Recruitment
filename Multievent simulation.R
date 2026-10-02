################################################################################
#          Simulation and Model Fitting for PCFG multi-event model             #
# Based on code from O. Gimenez https://github.com/oliviergimenez/banana-book/ #
################################################################################


#====================================
# Parameter and Population Settings
#====================================

# Population and Time -----------------
N       <- 300
time    <- 10
nstates <- 4 # 1-calf, 2-PCFG, 3-ENP, 4-dead/transient

# Parameters -------------------------------
# Initial state probabilities
initial_prob <- c(0.20, 0.60, 0.20) # Can't start in dead/transient 

# Transition probabilities
ps <- matrix(c(0.0, 0.7, 0.0, 0.3,
               0.0, 0.9, 0.0, 0.1,
               0.0, 0.2, 0.0, 0.8,
               0.0, 0.0, 0.0, 1.0),
             nrow = nstates, byrow = TRUE)

# Capture probabilities
po <- c(0.7, 0.80, 0.65, 0.0)

#===================================
# Simulation
#===================================
# packages
library(LaplacesDemon)

# Set-up data sets
f  <- sample(1:9,N,replace=TRUE)               # first observation in data set
Z  <- data.frame(matrix(NA,nrow=N,ncol=time))  # true underlying state
Y  <- data.frame(matrix(0,nrow=N,ncol=time))   # observations


for (i in 1:N){
  # Define latent state at first capture
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
# piP prob. of being in initial state PCFG
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

unascert <- nimbleCode({
  
  # Priors
  psiCP ~ dunif(0, 1)   # Prob transition of calf to PCFG (e.g., survival)
  psiPP ~ dunif(0, 1)   # Prob transition of PCFG to PCFG (e.g., survival)
  psiEP ~ dunif(0, 1)   # Prob transition of ENP to PCFG (e.g., permanent emigration)
  pC ~ dunif(0, 1)      # Prob observing PCFG calf
  pP ~ dunif(0, 1)      # Prob observing PCFG
  pE ~ dunif(0, 1)      # Prob observing ENP
  piC ~ dgamma(1, 1)    # Initial state calculation for calf
  piP ~ dgamma(1, 1)    # Initial state calculation for PCFG
  piE ~ dgamma(1, 1)    # Initial state calculation for ENP
  
  # vector of initial stats probs
  delta[1] <- piC / (piC + piP + piE) # prob. of being in initial state calf
  delta[2] <- piP / (piC + piP + piE) # prob. of being in initial state PCFG
  delta[3] <- piE / (piC + piP + piE) # prob of being in initial state ENP
  delta[4] <- 0                       # prob of being initial state dead
 
  
  # Transition probability matrix
  # probabilities of state z(t+1) given z(t)
  gamma[1,1:4] <- c(0, psiCP, 0, (1 - psiCP))
  gamma[2,1:4] <- c(0, psiPP, 0, (1 - psiPP))
  gamma[3,1:4] <- c(0, psiEP, 0, (1 - psiEP))
  gamma[4,1:4] <- c(0, 0, 0, 1)

  # Observation probability matrix
  # probabilities of y(t) given z(t)
  omega[1,1:5,1] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,1] <- c((1 - pP), 0, pP, 0, 0)
  omega[3,1:5,1] <- c((1 - pE), 0,  0, 0, pE)            
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
    # Define latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:4])
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
                omega.init = omega.init)


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
                                  piC = rgamma(1, 1, 1),
                                  piP = rgamma(1, 1, 1),
                                  piE = rgamma(1, 1, 1),
                                  z = zinit)}  

# MCMC settings
parameters.to.save <- c("psiCP",
                        "psiPP", 
                        "psiEP",
                        "pC",
                        "pP", 
                        "pE",
                        "delta")

n.iter   <- 500000
n.burnin <- 100000
n.chains <- 3
n.thin   <- 10


# Run the models
mcmc.unascert <- nimbleMCMC(code = unascert, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            thin = n.thin,
                            nburnin = n.burnin, 
                            nchains = n.chains)

MCMCsummary(mcmc.unascert, round = 2)

MCMCplot(mcmc.unascert)

MCMCtrace(object = mcmc.unascert,
          params = c("psiCP",
                     "psiPP", 
                     "psiEP",
                     "pC",
                     "pP", 
                     "pE",
                     "delta"))
