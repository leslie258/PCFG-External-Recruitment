################################################################################
#         Simulation and Model fitting for the JS PCFG multievent model        #
#           Initial State Probability a Function of Minimum Tenure             #
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
# Initial state probabilities
# Prob of being PCFG function of minimum tenure
# Prob of being ENP is 1 - pPCFG
# Prob of calf function of expected fecundity
init.C <- 0.24     # 0.24/1.24 = 0.20


# Probability of entry
Alpha <- 0.06

# Transition probabilities
ps <- matrix(c(0.0, 0.7, 0.0, 0.3, 0,
               0.0, 0.9, 0.0, 0.1, 0, 
               0.0, 0.2, 0.0, 0.8, 0,
               0.0, 0.0, 0.0, 1.0, 0),
               nrow = nstates-1, byrow = TRUE)


# Capture probabilities
po <- c(0.7, 0.8, 0.65, 0.0, 0.0)

#===================================
# Simulation
#===================================
# packages
library(LaplacesDemon)

# Set-up data sets
MT <- sample(1:182, N, replace=T)              # Minimum tenure
Z  <- data.frame(matrix(NA,nrow=N,ncol=time))  # true underlying state
Y  <- data.frame(matrix(0,nrow=N,ncol=time))   # observation

init.P <- numeric(0)
for (i in 1:N){
  # Define latent state at first capture
  Z[i,1] <- 5 # All individuals start as not yet entered
  for (t in 2:time){
    if(Z[i,t-1]==5){
      init.P[i] <- invlogit(-5 + 0.05*MT[i])
      init.E <- 1 - init.P[i]
      initial_prob <- c(init.C, init.P[i], init.E)/ (init.C + init.P[i] + init.E)
      ps.init <- c(initial_prob*Alpha, 0, (1-Alpha))
      Z[i,t] <- rcat(1,ps.init)
    }else{
    # State process: draw S(t) given S(t-1)
      Z[i,t] <- rcat(1, ps[Z[i,t-1], 1:5])}
    
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
# beta_0 intercept of prob of transition from ENP to PCFG
# beta_1 slope of prob of transition from ENP to PCFG
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

unascertJS.R <- nimbleCode({
  
  # Priors
  psiCP ~ dunif(0, 1)   # Prob transition of calf to PCFG (e.g., survival)
  psiPP ~ dunif(0, 1)   # Prob transition of PCFG to PCFG (e.g., survival)
  psiEP ~ dunif(0, 1)   # Prob transition of ENP to PCFG (e.g., permanent emigration)
  pC ~ dunif(0, 1)      # Prob observing PCFG calf
  pP ~ dunif(0, 1)      # Prob observing PCFG
  pE ~ dunif(0, 1)      # Prob observing ENP
  piC ~ dbeta(1, 3)     # Initial state calculation for calf
  beta_0 ~ dnorm(0,1)   # Intercept for prob PCFG model
  beta_1 ~ dnorm(0,1)   # Slope for prob PCFG model
  alpha ~ dunif(0, 1)   # Probability of entry
  
  
  # Observation matrix
  # probabilities of y(t) given z(t)
  omega[1,1:5,1] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,1] <- c((1 - pP), 0, pP, 0, 0)
  omega[3,1:5,1] <- c((1 - pE), 0,  0, 0, pE)            
  omega[4,1:5,1] <- c(1, 0, 0, 0, 0)
  omega[5,1:5,1] <- c(1, 0, 0, 0, 0)
  omega[1,1:5,2] <- c((1 - pC), pC, 0, 0, 0)
  omega[2,1:5,2] <- c((1 - pP),  0, 0, 0, pP)
  omega[3,1:5,2] <- c((1 - pE),  0, 0, 0, pE)            
  omega[4,1:5,2] <- c(1, 0, 0, 0, 0)
  omega[5,1:5,2] <- c(1, 0, 0, 0, 0)

  
  # Likelihood 
  for (i in 1:M){
    # All individuals start as having not entered the population
    z[i, 1] <- 5
    
    # Define initial state probabilities
    piP[i] <- ilogit(beta_0 + beta_1*MT[i])
    piE[i] <- 1 - piP[i]
    delta[1,i] <- piC/(piP[i] + piE[i] + piC)
    delta[2,i] <- piE[i]/(piP[i] + piE[i] + piC)
    delta[3,i] <- piP[i]/(piP[i] + piE[i] + piC)
    
    # Tranistion probability matrix
    # probabilities of state z(t+1) given z(t)
    gamma[1,1:5,i] <- c(0, psiCP, 0, (1 - psiCP), 0)
    gamma[2,1:5,i] <- c(0, psiPP, 0, (1 - psiPP), 0)
    gamma[3,1:5,i] <- c(0, psiEP, 0, (1 - psiEP), 0)
    gamma[4,1:5,i] <- c(0, 0, 0, 1, 0)
    gamma[5,1:5,i] <- c(alpha*delta[1,i], alpha*delta[2,i], alpha*delta[3,i], 0, (1-alpha))
    
    for (t in 2:K){
      # Determine timestep individual first seen
      first[i,t] <- 2 - equals(min(z[i,1:t]), 1) 
      
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:5,i])
      
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
my.data <- list(y = yy,
                MT = MT)


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
                                  piC = rbeta(1, 1, 3),
                                  beta_0 = rnorm(1, 0, 1),
                                  beta_1 = rnorm(1, 0, 1),
                                  alpha = runif(1, 0, 1),
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
                        "pE",
                        "alpha",
                        "N",
                        "Be")

n.iter   <- 500000
n.burnin <- 100000
n.chains <- 3
n.thin   <- 10


# Run the models
mcmc.unascertJSR <- nimbleMCMC(code = unascertJS.R, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            thin = n.thin,
                            nburnin = n.burnin, 
                            nchains = n.chains)

MCMCsummary(mcmc.unascertJSR, round = 2)

MCMCplot(mcmc.unascertJSR)

MCMCtrace(mcmc.unascertJSR)
          
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
MCMCsummary(mcmc.unascertJSR, round = 2, params = c('N[2]','N[3]','N[4]', 'N[5]', 'N[6]', 'N[7]', 'N[8]', 'N[9]', 'N[10]','N[11]'), 
            ISB = FALSE, exact = TRUE)

#=============================
# Entered from ENP
#=============================
enterENP <- function(x){
  sum(x==3)
}

apply(Z, 2, enterENP)
MCMCsummary(mcmc.unascertJSR, round = 2, 
            params = c('Be[2]','Be[3]','Be[4]', 'Be[5]', 'Be[6]', 'Be[7]', 'Be[8]', 'Be[9]', 'Be[10]','Be[11]'), 
            ISB = FALSE, exact = TRUE)

#==============================
# Model Parameters
#==============================

MCMCsummary(mcmc.unascertJSR, round = 2, 
            params = c("alpha","piC","pC","pP","pE","psiCP","psiPP","psiEP","beta_0","beta_1"), 
            ISB = FALSE, exact = TRUE)

# Alpha: 0.06
# pC: 0.7
# pP: 0.8
# pE: 0.65
# piC: 0.24
# psiCP: 0.7
# psiPP: 0.9
# psiEP: 0.2
# beta_0: -5
# beta_1: 0.05

