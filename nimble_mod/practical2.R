
library(nimble)
library(ggmcmc)
library(coda)

##############################################
##### Building  MCMC default algorithm #######
##############################################

# simulate some data

set.seed(24) 
M <- 10 
J <- 2 
y <- matrix(NA, nrow = M, ncol = J) 
psi <- 0.8 
p <- 0.5 
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

### declare NIMBLE/JAGS/BUGS code
occNimblecode <- nimbleCode({
  for (i in 1:nsite) {
    z[i] ~ dbern(psi) # True occupancy status
    for (j in 1:nvisit) {
      y[i,j] ~ dbern(z[i] * p) # Observed data
    }
  }
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  sum.z <- sum(z[1:nsite])
})


### Option 1 nimbleMCMC - all in one solution

nimbleMCMC_samples <- nimbleMCMC(code = occNimblecode, 
                                 constants = list(nvisit=J,nsite=M), 
                                 data = list(y=y), 
                                 inits = list(z = apply(y, 1, max),psi = 0.7, p = 0.33),
                                 nburnin = 1000, niter = 10000,monitors = c('psi','p'))

# returns samples as a matrix

par(mfrow=c(1,2))
plot(nimbleMCMC_samples[ , 'psi'], type = 'l', xlab = 'iteration',  
     ylab = expression(psi))
plot(nimbleMCMC_samples[ , 'p'], type = 'l', xlab = 'iteration',  
     ylab = expression(p))


# if I want to visualise mcmc samples with other packages
# specifically for mcmc
# I need ti fist declare nimbleMCMC_samples
# as an mcmc object using coda:as.mcmc
coda.samples <- as.mcmc(nimbleMCMC_samples)

coda.samples %>% ggs() %>% ggs_traceplot()



# Option 2 runMCMC (more steps but more flexible)


# step 1 create an operational from from NIMBLE/JAGS/BUGS code

occModel <- nimbleModel(code = occNimblecode, constants = list(nvisit=J,nsite=M),
                         data = list(y=y), inits = list(z = apply(y, 1, max),psi = 0.7, p = 0.33))


# step 2 build an MCMC object using buildMCMC(). we can add some customization here

occMCMC <- buildMCMC(occModel, monitors = c('psi','p'), thin = 3, useConjugacy = FALSE)

# useConjugacy controls whether conjugate samplers are assigned when possible

# you can now run the mcmc (Uncompiled MCMC will be really slow)

runMCMC(occMCMC, niter = 5) # not the best option


# step 3 before compiling the MCMC object we need to compile the model first

CoccModel <- compileNimble(occModel) # need to compile the model before I can compile the mcmc

# now the MCMC (project = NIMBLE model already associated with a project)

CoccMCMC <- compileNimble(occMCMC, project = occModel)

# and now we can use either $run or runMCMC() on the compiled model object.


runMCMC_samples <- runMCMC(CoccMCMC, nburnin = 1000, niter = 10000,nchains = 3,samplesAsCodaMCMC = T)

# last argument returns samples as mcmc object!

p1<- runMCMC_samples %>% ggs() %>% ggs_traceplot(greek = TRUE)
p2<- runMCMC_samples %>% ggs() %>% ggs_density(greek = TRUE)
gridExtra::grid.arrange(p1,p2,ncol=2)

# this option is the same, but I cannot assign the results to an object and
# I don't have the same number of srguments to custimise the run call

CoccMCMC$run(10000)
run_samples <- as.matrix(CoccMCMC$mvSamples)

plot(runMCMC_samples[ , 'psi'], type = 'l', xlab = 'iteration',  ylab = expression(psi))
plot(runMCMC_samples[ , 'p'], type = 'l', xlab = 'iteration',  ylab = expression(p))



#############################################################################################
######## Defining alternative models sensitivity analysis and a tiny simulation study #######
#############################################################################################

##### Alternative models

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


set.seed(123)
occ.data <- sim.occ(mean.psi = 0.8,p = 0.4,beta = c(0.5,-1.5,1.5),J = 3,M = 10)


# define de model

occNimblecovariates <- nimbleCode({
  
  # Priors
  # Intercept
  mean.psi ~ dunif(0, 1) # Detection intercept on prob. scale
  beta0 <- logit(mean.psi) # Detection intercept
  p ~ dunif(0,1)
  
  # alternative model priors for beta
  
  for (k in 1:ncov){
    
    if(beta.prior == "logistic"){
      beta[k] ~ dlogis(0,1)
    } else if(beta.prior == "t-prior"){
      beta[k] ~ dt(0,1.566,7.763)
    } else if(beta.prior == "w-normal"){
      beta[k] ~ dnorm(0, pow(2.25,-2))
    } else{
      beta[k] ~ dnorm(0, pow(1.4,-2))
    }
  }

  # Likelihood
  for (i in 1:nsite) {
    z[i] ~ dbern(psi[i]) 
    logit(psi[i]) <- beta0 + inprod(beta[1:ncov], x[i, 1:ncov])
    # beta[1] * x[i,1] + beta[2] * x[i,2] + ...
    # beta[1:ncov] %*% x[i, 1:ncov]
    
    for (j in 1:nvisit) {
      y[i,j] ~ dbern(z[i] * p) 
    }
  }
  
})


# declared the argument within the if-then-else statement to define alternative models, e.g. lets compare two priors

beta.prior <- "w-normal"
model.normalprior <- nimbleModel(code = occNimblecovariates,data = list(y=occ.data[[1]],x=occ.data[[2]]),
                           constants = list(nvisit=3,nsite=10,ncov=3),
                           inits = list(mean.psi = runif(1),
                                        z = apply(occ.data[[1]], 1, max),
                                        beta =rnorm(3), p = runif(1)))

beta.prior <- "t-prior"
model.tprior <- nimbleModel(code = occNimblecovariates,data = list(y=occ.data[[1]],x=occ.data[[2]]),
                           constants = list(nvisit=3,nsite=10,ncov=3),
                           inits = list(mean.psi = runif(1),
                                        z = apply(occ.data[[1]], 1, max),
                                        beta =rnorm(3), p = runif(1)))


model.tprior$getCode()
model.normalprior$getCode()

mcmc1 <- buildMCMC(model.normalprior, monitors = c('beta'))
mcmc2 <- buildMCMC(model.tprior, monitors = c('beta'))

Cmodel1 <- compileNimble(model.normalprior)
Cmodel2 <- compileNimble(model.tprior)

Cmcmc1<- compileNimble(mcmc1,project = model.normalprior)
Cmcmc2<- compileNimble(mcmc2,project = model.tprior)
  
  
runMCMC_samples1 <- runMCMC(Cmcmc1, nburnin = 1000, niter = 10000,nchains = 3,samplesAsCodaMCMC = T)
runMCMC_samples2 <- runMCMC(Cmcmc2, nburnin = 1000, niter = 10000,nchains = 3,samplesAsCodaMCMC = T)

p1<- runMCMC_samples1 %>% ggs() %>% ggs_traceplot(greek = TRUE)
p2<- runMCMC_samples2 %>% ggs() %>% ggs_traceplot(greek = TRUE)
gridExtra::grid.arrange(p1,p2,ncol=2)


####### Simulation study ########

# compare the MSE with varying values of p


pdet<-c(0.1,0.25,0.5,0.75)
MSE<-matrix(NA,ncol=4,nrow=20)

for(d in 1:4){
  for(r in 1:20){
    occ.data <- sim.occ(mean.psi = 0.8,p = pdet[d],beta = c(0.5,-1.5,1.5),J = 3,M = 10)
    Cmodel1$y<- occ.data[[1]]
    Cmodel1$x<- occ.data[[2]]
    Cmodel1$z<- apply(occ.data[[1]], 1, max)
    samples <- runMCMC(Cmcmc1, nburnin = 1000, niter = 10000,samplesAsCodaMCMC = T)
    MSE[r,d]<- apply(samples,1, function(x){(x-c(0.5,-1.5,1.5))%*%diag(3)%*%(x-c(0.5,-1.5,1.5))}) %>% mean()
    
  }
}

ggplot(data=reshape::melt(MSE),aes(x=as.factor(X2),y=value,fill=as.factor(X2)))+
  geom_boxplot()+
  labs(y='MSE',x='p')+
  scale_x_discrete(labels=c("1" = "0.1", "2" = "0.25", "3" = "0.5", "4" = "0.75"))+
  theme(legend.position = 0)

###################################################
############# Customizing a MCMC ##################
###################################################


occ.data <- sim.occ(mean.psi = 0.8,p = 0.5,beta = c(0.5,-1.5,1.5),J = 3,M = 100)

## constants, data, and initial values
constants <- list(nvisit=3,nsite=10,ncov=3)

data <- list(
  y=occ.data[[1]],
  x=occ.data[[2]]
)

inits <- list(mean.psi = runif(1),
              z = apply(occ.data[[1]], 1, max),
              beta =rnorm(3),
              p = runif(1))

beta.prior <- "w-normal"
model.normalprior <- nimbleModel(code = occNimblecovariates,data = data,
                                 constants = constants,
                                 inits = inits)

mcmcOcc <- buildMCMC(model.normalprior)

MCMCconfiguration <- configureMCMC(model.normalprior,monitors = 'beta',thin = 3)
# use this if you want to customise the thinning for each parameter
# otherwise just use the monitor parameter in the previous functions

# you can include two sets of monitors, each with different thinning intervals.
MCMCconfiguration$addMonitors2('p') 
MCMCconfiguration$setThin2(1)


MCMCconfiguration$printSamplers()

# beta[1] through beta[3] have each been assigned adaptive random walk Metropolis-Hastings samplers.

MCMCconfiguration$monitors;MCMCconfiguration$monitors2

# replace the univariate samplers with a block sampler
# RW block sampler performs a simultaneous update of one or more model nodes, using
# the Metropolis-Hastings algorithm with a multivariate normal proposal distribution
# Automated Factor Slice Sampler 'AF_slice' performs better in practice

# step 1 remove old sampler
MCMCconfiguration$removeSamplers('beta', print = FALSE)
# step 2 add new sampler 
MCMCconfiguration$addSampler(target = 'beta[1:3]', type = 'RW_block')

# Build the customized MCMC
customMCMC <- buildMCMC(MCMCconfiguration)

# Compile the model and MCMC

Custom_model_compiled <- compileNimble(model.normalprior, customMCMC)

Custom_model_compiled$customMCMC$run(10000,time = TRUE)
# a vector of the total time spent in each sampler, measured in seconds.
Custom_model_compiled$customMCMC$getTimes()


MCMCsamples <- as.matrix(Custom_model_compiled$customMCMC$mvSamples)
pairs(MCMCsamples, pch = '.')

plot(as.matrix(Custom_model_compiled$customMCMC$mvSamples2), 
     type = 'l', xlab = 'iteration',  ylab = expression(p))

###########################################################
################# End of practical 2 ######################
###########################################################



