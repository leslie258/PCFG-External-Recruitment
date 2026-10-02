################################################################################
#         Simulation and Model fitting for the JS PCFG multievent model        #
# Based on code from O. Gimenez https://github.com/oliviergimenez/banana-book/ #
################################################################################


#====================================
# Parameter and Population Settings
#====================================

# Population and Time -----------------
N       <- 600
time    <- 11
nstates <- 5 # 1-calf, 2-PCFG, 3-ENP, 4-dead/transient, 5-not yet entered


# Parameters -------------------------------
# Probability of entry (a)
Alpha <- 0.06

# Transition probabilities
ps <- matrix(c(0.0, 0.7, 0.0, 0.3, 0,
               0.0, 0.9, 0.0, 0.1, 0, 
               0.0, 0.2, 0.0, 0.8, 0,
               0.0, 0.0, 0.0, 1.0, 0,
               Alpha*0.2, Alpha*0.6, Alpha*0.2, 0, (1-Alpha)),
               nrow = nstates, byrow = TRUE)


# Capture probabilities
po <- c(0.7, 0.8, 0.65, 0.0, 0.0)

#===================================
# Simulation
#===================================
# packages
library(LaplacesDemon)

# Set-up data sets
Z  <- data.frame(matrix(NA,nrow=N,ncol=time))  # true underlying state
Y  <- data.frame(matrix(0,nrow=N,ncol=time))   # observation


for (i in 1:N){
  # Define latent state at first capture
  Z[i,1] <- 5 # All individuals start as not yet entered
  for (t in 2:time){
    # State process: draw S(t) given S(t-1)
    Z[i,t] <- rcat(1, ps[Z[i,t-1], 1:5])
    # Observation process: draw O(t) given S(t)
    Obs   <- rbinom(1, 1, po[Z[i,t]])
    if(Z[i,t]==5){
      Zt <- 0
    }else{
      f  <- which(Z[i,1:t] == min(Z[i,1:t]))[1]
      if(f==t){
        Zt <- ifelse(Z[i,f]==1, 1, 3)
      }else{
      Zt <- ifelse(Y[i,f]==1, 2, 3)}
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

#====================================
# Nimble Implementation
#====================================

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
# alpha prob of entry
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
# 5 = not yet entered
# -------------------------------------------------

unascertJS <- nimbleCode({
  
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
  alpha ~ dunif(0, 1)   # Probability of entry
  
  # vector of initial stats probs
  delta[1] <- piC / (piC + piP + piE) # prob. of being in initial state calf
  delta[2] <- piP / (piC + piP + piE) # prob. of being in initial state PCFG
  delta[3] <- piE / (piC + piP + piE) # prob of being in initial state ENP
  delta[4] <- 0                       # prob of being in initial state dead
  delta[5] <- 0                       # prob of being in initial state not entered
 
  
  # Transition probability matrix
  # probabilities of state z(t+1) given z(t)
  gamma[1,1:5] <- c(0, psiCP, 0, (1 - psiCP), 0)
  gamma[2,1:5] <- c(0, psiPP, 0, (1 - psiPP), 0)
  gamma[3,1:5] <- c(0, psiEP, 0, (1 - psiEP), 0)
  gamma[4,1:5] <- c(0, 0, 0, 1, 0)
  gamma[5,1:5] <- c(alpha*delta[1], alpha*delta[2], alpha*delta[3], 0, (1-alpha))

  # Observation probability matrix
  # probabilities of y(t) given z(t)
  omega[1,1:5,1] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,1] <- c((1 - pP), 0, pP, 0, 0)
  omega[3,1:5,1] <- c((1 - pE), 0,  0, 0, pE)            
  omega[4,1:5,1] <- c(1, 0, 0, 0, 0)
  omega[5,1:5,1] <- c(1, 0, 0, 0, 0)
  omega[1,1:5,2] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,2] <- c((1 - pP), 0,  0, 0, pP)
  omega[3,1:5,2] <- c((1 - pE), 0,  0, 0, pE)            
  omega[4,1:5,2] <- c(1, 0, 0, 0, 0)
  omega[5,1:5,2] <- c(1, 0, 0, 0, 0)

  
  # Likelihood 
  for (i in 1:M){
    # All individuals start as having not entered the population
    z[i, 1] <- 5
    #y[i,first[i]] <- omega.init[z[i,first[i]]]
    for (t in 2:K){
      # Determine timestep individual first seen
      first[i,t] <- 2 - equals(min(z[i,1:t]), 1) 
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:5])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:5,first[i,t]])
    }
  }
  
  # Derived population parameters
  for (i in 1:M){
    # Alive at time t-1
    for (t in 1:K){
      alC[i,t] <- equals(z[i,t], 1)
      alP[i,t] <- equals(z[i,t], 2)
      alE[i,t] <- equals(z[i,t], 3)
    } #t 
    
    # Entered the population at time t
    for (t in 2:K){
      dC[i,t] <- equals((z[i,t-1]-alC[i,t]),4)*alC[i,t]
      dP[i,t] <- equals((z[i,t-1]-alP[i,t]),4)*alP[i,t]
      dE[i,t] <- equals((z[i,t-1]-alE[i,t]),4)*alE[i,t]
    }
    
   # Whether individual was alive in any of the states during the time series  
    aliveC[i] <- sum(alC[i,1:K]) #alive as a calf
    aliveP[i] <- sum(alP[i,1:K]) #alive as PCFG
    aliveE[i] <- sum(alE[i,1:K]) #alive as ENP
    
  } # i

  # Abundance
  for (t in 1:K){
    Nc[t] <- sum(alC[1:M,t])  # Number of calves at time t
    Np[t] <- sum(alP[1:M,t])  # Number of PCFG at time t
    Ne[t] <- sum(alE[1:M,t])  # Number of ENP at time t
    N[t]  <- sum(alC[1:M,t]) + sum(alP[1:M,t]) + sum(alE[1:M,t])} # Actual population size
    
  for(t in 2:K){
    Bc[t] <- sum(dC[1:M,t])  # Number of calf entries
    Bp[t] <- sum(dP[1:M,t]) # Number of PCFG entries
    Be[t] <- sum(dE[1:M,t])  # Number of ENP entries
    B[t]  <- sum(dC[1:M,t]) + sum(dP[1:M,t]) + sum(dE[1:M,t]) # Number of entries
  } #t
  
  # Superpopulation
  for (i in 1:M){
    alive[i] <- aliveC[i] + aliveP[i] + aliveE[i]
    w[i] <- 1 - equals(alive[i],0)
  } #i
  
  Nsuper <- sum(w[1:M]) # Superpopulation size

  }) # end Nimble code
  


# constants
get.first <- function(x){
  if(sum(x)==0){
    1
  }else{
    min(which(x != 0))}
}
first <- apply(y, 1, get.first)
my.constants <- list(K = ncol(y), 
                     M = nrow(y))

# data
yy <- y + 1
yy[yy==4] <- 5
my.data <- list(y = yy)


# Initial values for Z
zinit <- y

for(i in 1:nrow(y)) {
  for (j in 1:ncol(y)){
    if(j == first[i] & y[i,first[i]]==3) {zinit[i,j] <- sample(c(2,3),1)}
    if(j > first[i]) {zinit[i,j] <- 2}
    if(first[i] == 1 | j < first[i]) {zinit[i,j] <- 5}
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
                                  alpha = runif(1, 0, 1),
                                  z = zinit)}  

# MCMC settings
parameters.to.save <- c("psiCP",
                        "psiPP", 
                        "psiEP",
                        "pC",
                        "pP", 
                        "pE",
                        "alpha",
                        "delta",
                        "N",
                        "Be")

n.iter   <- 500000
n.burnin <- 100000
n.chains <- 3
n.thin   <- 10


# Run the models
mcmc.unascertJS <- nimbleMCMC(code = unascertJS, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            thin = n.thin,
                            nburnin = n.burnin, 
                            nchains = n.chains)


MCMCplot(mcmc.unascertJS)

MCMCtrace(object = mcmc.unascertJS)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Check results against simulated data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#===========================
# Abundance
#===========================

abundance <- function(x){
  n <- x[x != 5 & x != 4]
  length(n)
  
}

apply(Z, 2, abundance)
MCMCsummary(mcmc.unascertJS, round = 2, params = c('N[2]','N[3]','N[4]', 'N[5]', 'N[6]', 'N[7]', 'N[8]', 'N[9]', 'N[10]','N[11]'), 
            ISB = FALSE, exact = TRUE)


#=============================
# Entered from ENP
#=============================
enterENP <- function(x){
  sum(x==3)
}

apply(Z, 2, enterENP)
MCMCsummary(mcmc.unascertJS, round = 2, 
            params = c('Be[2]','Be[3]','Be[4]', 'Be[5]', 'Be[6]', 'Be[7]', 'Be[8]', 'Be[9]', 'Be[10]','Be[11]'), 
            ISB = FALSE, exact = TRUE)

#==============================
# Model Parameters
#==============================

MCMCsummary(mcmc.unascertJS, round = 2, 
            params = c("alpha","pC","pP","pE","psiCP","psiPP","psiEP","delta[1]","delta[2]","delta[3]"), 
            ISB = FALSE, exact = TRUE)

# Alpha: 0.06
# pC: 0.7
# pP: 0.8
# pE: 0.65
# psiCP: 0.7
# psiPP: 0.9
# psiEP: 0.2
# delta[1]: 0.2
# delta[2]: 0.6
# delta[3]: 0.2


