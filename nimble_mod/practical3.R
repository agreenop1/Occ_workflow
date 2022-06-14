################################################################
################ Posterior Predictive Checks ###################
################################################################

sim.occ <- function(mean.psi,p,beta,J,M){
  
  nsites <- M
  beta0 <- qlogis(mean.psi)
  ncov<-length(beta)
  X <- matrix(rnorm(M*ncov), nrow = M, ncol = ncov)
  psi.occ <- plogis(beta0 +X%*%beta)
  z <- rbinom(n = nsites, size = 1, prob = psi.occ)
  ysim <- matrix(NA, nrow = nsites, ncol = J)
  for(j in 1:J){
    ysim[,j] <- rbinom(n = nsites, size = 1, prob = z*p)
  }
  return(list(ysim,X))
}


occ.data <- sim.occ(mean.psi = 0.8,p = 0.5,beta = c(0.5,-1.5,1.5),J = 3,M = 500)

## constants, data, and initial values
constants <- list(V=3,nsite=500,ncov=3,zeros = rep(0, 3), omega = 0.0001 * diag(3))

data <- list(
  y=apply(occ.data[[1]], 1, sum),
  x=occ.data[[2]]
)

inits <- list(mean.psi = runif(1),
              z = apply(occ.data[[1]], 1, max),
              beta =rnorm(3),
              p = runif(1))


# The code 
Modelcode <- nimbleCode({
  
  # Priors
  # Intercept
  mean.psi ~ dunif(0, 1) 
  beta0 <- logit(mean.psi) 
  p ~ dunif(0,1)
  
  beta[1:ncov] ~ dmnorm(zeros[1:ncov], omega[1:ncov, 1:ncov])
  linpred[1:nsite] <- (x[1:nsite, 1:ncov] %*% beta[1:ncov])[1:nsite,1]
  
  # Likelihood
  for (i in 1:nsite) {
    z[i] ~ dbern(psi[i]) 
    logit(psi[i]) <- beta0 + linpred[i]
    y[i] ~ dbin(p * z[i], V) 
    # GOF assessment Freeman-Tukey discrepancy
     FT[i] <- (sqrt(y[i]) - sqrt(p*z[i]*V))^2 
  }
  
  Tobs<-sum(FT[1:nsite])
  
})

# Operational model

OccModel <- nimbleModel(Modelcode, constants = constants,
                        inits = inits,data=data)



## Ensure we have the nodes needed to simulate new datasets
dataNodes <- OccModel$getNodeNames(dataOnly = TRUE)
parentNodes <- OccModel$getParents(dataNodes, stochOnly = TRUE)
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- OccModel$getDependencies(parentNodes, self = FALSE)



OccModel.compiled <- compileNimble(OccModel)
mcmc    <- buildMCMC(OccModel,monitors = c(parentNodes,"Tobs"))
cmcmc   <- compileNimble(mcmc, project = OccModel)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)


nSamp <- nrow(samples)
ppSamples <- numeric(nSamp)
vars <- cmcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix

set.seed(1)
system.time({
for(i in 1:nSamp){
  values(OccModel.compiled, vars) <- samples[i, ]  # put values from a vector into a set of model nodes
  OccModel.compiled$simulate(simNodes, includeData = TRUE)
  ppSamples[i] <- values(OccModel.compiled,"Tobs")
  
}})

#  user  system elapsed 
# 657.85    0.00  657.89 


# plot(samples[,"Tobs"],ppSamples)
# abline(0,1,col=2,lwd=3)

## too slow! lets create a nimble function instead

ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)   
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, 'FT')
    }
    returnType(double(2))       
    return(ppSamples)
  })

ppSampler <- ppSamplerNF(OccModel, mcmc)
cppSampler <- compileNimble(ppSampler, project = OccModel)
identical(colnames(samples), OccModel$expandNodeNames(mcmc$mvSamples$getVarNames()))

samples<- samples[, OccModel$expandNodeNames(mcmc$mvSamples$getVarNames())]
identical(colnames(samples), OccModel$expandNodeNames(mcmc$mvSamples$getVarNames()))

set.seed(1)
system.time(ppSamples_via_nf <- cppSampler$run(samples))

# user  system elapsed 
# 0.31    0.03    0.35 

Tsim <- apply(ppSamples_via_nf,1,sum)


plot(samples[,"Tobs"],apply(ppSamples_via_nf,1,sum),ylab='Simulated data',xlab='Observed data')
abline(0,1,col=2,lwd=3)
bpval <- round(mean(apply(ppSamples_via_nf,1,sum)
> samples[,"Tobs"]),digits=2)
unikn::mark(labels = paste("Bayesian p-val =", bpval), x = 45, y = 110,cex=1.3)



df<- data.frame(TF = c(samples[,"Tobs"],Tsim),
                values = rep(c('fitted','simulated'),each=nSamp))

library(ggplot2)
ggplot(df, aes(x=TF, fill=values)) +
  geom_density(alpha=0.4)



#####################################################################
################## user-defined distributions #######################
#####################################################################

dZIB <- nimbleFunction(
  run = function(x = integer(), p = double(), size  = integer(),
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dbinom(x, size, p, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dbinom(x, size, p, log = FALSE))
    }
    ## when x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dbinom(0, size, p, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })


rZIB <- nimbleFunction(
  run = function(n = integer(), p = double(), size = integer(), zeroProb = double()) {
    returnType(integer())
    z <- rbinom(1, size=1, prob=zeroProb)
    y <- (1-z)*rbinom(1,size=size,prob=p)
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(y)
  })

ZIBcode <- nimbleCode({
  
  # Priors
  psi ~ dunif(0, 1) 
  p ~ dunif(0,1)
  phi <- 1 - psi
  
  for (i in 1:nsite){
    y[i] ~ dZIB(size = V,p = p, zeroProb = phi) 
  }
  
  # ZIB Likelihood to derive the Deviance
  
  for( i in 1:nsite){
    BinCoef[i] <- exp(lfactorial(V) - (lfactorial(y[i]) + lfactorial(V-y[i])))
    lik[i] <- BinCoef[i]*psi*(p^y[i])*(1-p)^(V-y[i]) + (y[i]==0)*(1-psi)
  }
  
  Deviance <- -2*log(prod(lik[1:nsite]))
  
})



set.seed(24) 
M <- 100 
J <- 4
y <- matrix(NA, nrow = M, ncol = J) 
psi <- 0.8 
p <- 0.5 
z <- rbinom(n = M, size = 1, prob = psi) 
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}
ysum<- apply(y,1,sum)

samples <- nimbleMCMC(code = ZIBcode, WAIC = T,summary = T,
                      constants = list(V=J,nsite=M), 
                      data = list(y=ysum), 
                      inits = list(psi = 0.7, p = 0.33),
                      nburnin = 1000, niter = 10000,monitors = c('psi','p','Deviance'))

samples$summary
samples$WAIC

par(mfrow=c(1,2))
plot(samples$samples[ , 'psi'], type = 'l', xlab = 'iteration',  
     ylab = expression(psi))
plot(samples$samples[ , 'p'], type = 'l', xlab = 'iteration',  
     ylab = expression(p))



# R Deviance of posterior mean and DIC calculation
lik.post<-numeric(M)
for(i in 1:M){
  psi <- samples$summary[3,1]
  p <-  samples$summary[2,1]
  lik.post[i] <- choose(J,ysum[i])*psi*(p^ysum[i])*(1-p)^(J-y[i]) + (y[i]==0)*(1-psi)
}

Post.dev <- -2*log(prod(lik.post))
pD <- samples$summary[1,1] - Post.dev
DIC <- samples$summary[1,1] + pD
DIC
pV <- (samples$summary[1,3]**2)/2
DICjags<- Post.dev + pV

# DIC = posterior mean of the deviance +  the deviance evaluated at the posterior mean of the model parameter(s)




 
####################################################################
################ Nimble Ecology Example ############################
####################################################################

library(nimbleEcology)


occupancy_code <- nimbleCode({
  psi ~ dunif(0,1)
  p ~ dunif (0,1)
  for(i in 1:nSites) {
    y[i, 1:nVisits] ~ dOcc_s(probOcc = psi, probDetect = p, len = nVisits)
  }
})

occ_samples<- nimbleMCMC(code = occupancy_code, summary = T,
           data = list(y = y),
           constants = list(nSites = M, nVisits = J),
           inits = list(psi = 0.7, p = 0.15),
           nburnin = 1000, niter = 10000,monitors = c('psi','p'))

occ_samples$summary


########################################
############### The End ################
########################################




