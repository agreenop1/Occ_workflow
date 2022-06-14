############################################################################################
# 2. Multi-species occupancy model for risk quotient & environment with random species     #
# effect                                                                                   #
############################################################################################
rm(list = ls())
library(nimble)


# slurm input
i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
source('nimble_save_f.R')
set.seed(i) # 1:3 depending on chain

group = 'hoverflies_pol'
mod = 3
cat(mod,group,"\n")

# data formatted for jags


# output files
rdataFile=paste0("Model_outputs/",mod,"_",group,"_all_C.",i,"_run",".rdata")
save.dir=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",mod,"_",group)

jas_data <- readRDS("Model_data/data_hoverflies_all.499_1994.2016pol.rds")


# CHECK MCMC PARAMETERS !! # 
# MCMC settings
ni <- 100; nt <- 5 ; n.species="2"; daisy=T # can be all

# covariate formula
covs <- jas_data[[4]]# covariates

cat(names(covs),"\n")
######################################################################
# data in format from jags data prep

occup<- jas_data[[1]] #occ data
visit<- jas_data[[2]]# visit info
zobs<- jas_data[[3]]# init values
closure.period <- jas_data[[5]] #closure period

# check sample sizes in original data set
nsites <- length(unique(visit$site_5km)) # sites
nclosure <- length(unique(visit$TP)) # number of time points
nspecies <- ncol(occup[-1]) # num. species

# CHECK SPECIES NUMBER !
nspecies = if(n.species=="all"){nspecies}else{n.species}
zobs <- zobs[1:nspecies,,] # subset chosen for test if applic.


# remove id column
y <- occup[-1]

###################################################################
# covarite processing

# create array for covs
row = nrow(covs[[1]]); col = ncol(covs[[1]]); nmat = length(covs)
cov.array <- array(dim=c(row,col,nmat))

for(n in 1:nmat){
  cov.array[,,n]  <- as.matrix(covs[[n]])
}


##########################################################################
################# OCCUPANCY MODEL ########################################

# Bundle data and summarize data bundle
win.data <-    list(y = y[1:nspecies] , 
                    closure = visit$TP, site = visit$site_5km.n,
                    nclosure = nclosure, nsites = nsites, nspecies = nspecies,
                    SHORT = visit$SHORT, 
                    LONG = visit$LONG,
                    nobs = nrow(occup),
                    COVS=cov.array,
                    mod=mod)

# Specify model in BUGS language for vertical data format
modx  <- nimbleCode({
  ######################### PRIORS ##################################
  # ECOLOGICAL MODEL #
  # by species
  for (s in 1:nspecies){ 
    
    logitinit[s] ~ dnorm(mu.init,tau.init)  #  initial occ prior -heterogeneity between species
    logit(init.occ[s]) <- logitinit[s]
   
    
    logitgamma[s] ~ dnorm(mu.gamma,tau.gamma)  #  colonisation prior -heterogeneity between species
    logit(gamma[s]) <- logitgamma[s]
    
    alpha.phi[s] ~ dnorm(mu.alpha.phi,tau.alpha.phi) # intercept for species persistence model logit
    
    for (n in 1:ncov){  
      beta[n,s] ~ dnorm(mu.beta[n],tau.beta[n]) # covariate array prior
    }
    
  }
  
  # ecology model hyperpriors
  mu.alpha.phi ~ dnorm(0,0.01) # persistence intercept (expected value=0, so persistence=0.5)
  mu.gamma ~ dnorm(0,0.01) # colonisation intercept
  mu.init ~ dnorm(0,0.01) # init intercept
  
  tau.alpha.phi ~ T(dt(0,1,1),0,) # for SD truncated student t - approx half cauchy
  tau.gamma ~ T(dt(0,1,1),0,)
  tau.init ~ T(dt(0,1,1),0,)
  
  # covariate hyperpriors
  
  if(mod==1){
    ncov <- 6 # linear - separate spatial and temporal term
  } else if (mod==2){
    ncov <- 5 # quadratic term for temperature (no different spatial/tempor effects for temp or rq)
  } else {
    ncov <- 8 #  quadratic term separate spatial and temporal effects
  }
  
  for(n in 1:ncov){
    beta.mean[n] ~ dunif(0,1) 
    mu.beta[n] <- logit(beta.mean[n]) # mean
    tau.beta[n] ~ T(dt(0,1,1),0,)
  }
  
  
  # OBSERVATION MODEL
  # Observation model priors
  for(i in 1:nspecies){
    dtype1.p[i] ~ dnorm(mu.d1.p,tau.p1)
    dtype2.p[i] ~ dnorm(mu.d2.p,tau.p2)
    dtype3.p[i] ~ dnorm(mu.d3.p,tau.p3)
  }
  
  for (t in 1:nclosure) {
    alpha.p[t] ~ dnorm(0,tau.p)
  }
  # observation model hyperpriors 
  mu.d1.p ~ dnorm(-2,0.01) # mean species effect = 0.12 (probability of being recorded) on list length 1
  mu.d2.p ~ dnorm(0,0.01)
  mu.d3.p ~ dnorm(0,0.01)
  
  tau.p1 ~ T(dt(0,1,1),0,)
  tau.p2 ~ T(dt(0,1,1),0,)
  tau.p3 ~ T(dt(0,1,1),0,)
  
  tau.p ~ T(dt(0,1,1),0,) # for between year variance#
  
  ##################### END OF PRIORS #########################
  
  
  # Ecological model: Define state conditional on parameters
  for (s in 1:nspecies){ 
    for (i in 1:nsites){ # 5km site
      z[s,i,1] ~ dbern(init.occ[s])
      for (t in 2:nclosure){
        # persistence at i for t is function of pesticides
        logit(phi[s,i,t]) <- 
          
          if(mod==1){
            
          alpha.phi[s] + 
          
          beta[1,s]*COVS[i,t-1,1]+ 
          beta[2,s]*COVS[i,t-1,2]+ 
          beta[3,s]*COVS[i,t-1,3]+ 
          beta[4,s]*COVS[i,t-1,4]+ 
          beta[5,s]*COVS[i,t-1,5]+ 
          beta[6,s]*COVS[i,t-1,6]
          
          } else if (mod==2){
            
          alpha.phi[s] + 
          
          beta[1,s]*COVS[i,t-1,7]+ 
          beta[2,s]*pow(COVS[i,t-1,7],2)+ 
          beta[3,s]*COVS[i,t-1,3]+ 
          beta[4,s]*COVS[i,t-1,4]+ 
          beta[5,s]*COVS[i,t-1,8] 
        
          
          } else {
            
          alpha.phi[s] + 
          
          beta[1,s]*COVS[i,t-1,1]+ 
          beta[2,s]*pow(COVS[i,t-1,1],2)+ 
          beta[3,s]*COVS[i,t-1,2]+ 
          beta[4,s]*pow(COVS[i,t-1,2],2)+ 
          beta[5,s]*COVS[i,t-1,3]+ 
          beta[6,s]*COVS[i,t-1,4]+ 
          beta[7,s]*COVS[i,t-1,5]+ 
          beta[8,s]*COVS[i,t-1,6]
          
          }
        
        z[s,i,t] ~ dbern(z[s,i,t-1]*phi[s,i,t] + (1-z[s,i,t-1])*gamma[s])
      }
    }
  }
  
  # Observation model
  for (s in 1:nspecies){
    for (i in 1:nobs){
      logit(p[s,i]) <- alpha.p[closure[i]] + dtype1.p[s] + dtype2.p[s]*SHORT[i] +
        dtype3.p[s]*LONG[i]
      
      y[i,s] ~ dbern(z[s,site[i],closure[i]]*p[s,i])
    }
  }
})




# Initial values
zst <- zobs
inits <-list(z = zst)

# Parameters monitored
params <-  c("beta","mu.beta","tau.beta","z","mu.alpha.phi","mu.gamma","tau.beta",
             "deviance","alpha.p","phi",
             "dtype1.p","dtype2.p","dtype3.p","init.occ","alpha.phi","gamma","mu.init")  

# Call nimble from R

mod. <- nimbleModel(code = mod, name = "mod.",  constants = win.data,
                    inits = inits)

# build mcmc - in this case all the default samplers are used
modmcmc <- buildMCMC(mod.,monitors=params,thin=nt)

# compile nimble model
cmod <- compileNimble(mod.)

# compile mcmc algorithm
cmcmc <- compileNimble(modmcmc,project = mod.)

# run the algorithm 
cmcmc$run(niter=100,reset=F,resetMV = T)

# save the state of the model
stateList <- list(modelState = getModelState(cmod),
                  mcmcState = getMCMCstate(modcon, cmcmc))  


# daisy ID
id <- paste0("_ID_",1)

# remove samples
out <- list(samples=as.matrix(cmcmc$mvSamples),parameters=params,modState=stateList,iterations=ni,modData=win.data,inits=inits,
            chain=i,id=id)

# save daisy
if(daisy){
  save(out,file=paste0(save.dir,"_C.",i,id,".rdata"))
}

# save runtime
save(out,file=rdataFile)

