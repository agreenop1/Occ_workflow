############################################################################################
# Multi-species occupancy model for risk quotient & environment with random species effect #
############################################################################################
# individual RQ analysis - no environmental variables
rm(list = ls())

library(tictoc)

i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]


# data formatted for jags
jas_data <- readRDS("Model_data/data_spiders_all.499_1994.2010.rds")

# output files
rdataFile=paste0("Model_outputs/mod.spi.B.all_C.",i,"_run",".rdata")
save.dir="Model_outputs/mod.spi.B.all"

# CHECK MCMC PARAMETERS !! # 
# MCMC settings
ni <- 100; nt <- 3 ; nb <- 0 ; nc <- 1; na <- 100; n.species="all"; daisy=T # can be all
n.species=2

# covariate formula 
covars<- jas_data[[4]]# covariates 

covs<-list(covars$temp,covars$rq_sum_wide,covars$semi)
inters <-NULL

######################################################################
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
y[y==T] <- 1

# if closure period is 1 then earliest year = 1 etc.
# ensures data is compatible with jags
if(closure.period==1){visit$TP <- visit$TP-(min(years)-1) }


##########################################################################
################# OCCUPANCY MODEL ########################################
row = nrow(covars[[1]]); col = ncol(covars[[1]]); nmat = length(covs)
cov.array <- array(dim=c(row,col,nmat))

for(n in 1:nmat){
  cov.array[,,n]  <- as.matrix(covs[[n]])
}

formula <- paste0("beta[",1:nmat,",s]*COVS[i,t-1,",1:nmat,"]",collapse ="+")


if(!is.null(inters)){
  coms <- combn(inters,2)  
  
  for (n in 1:ncol(coms)){
    ints1 <- paste0("+beta[",n+nmat,",s]*COVS[i,t-1,",coms[1,n],"]*COVS[i,t-1,",coms[2,n],"]")
    formula <- paste0(formula,ints1) }

nmat <- ncol(coms)+nmat
  }

vars <- strsplit(formula,split="+",fixed=T)[[1]]

cat("Ecological model covariate structure",vars,sep="\n")


# Bundle data and summarize data bundle
  win.data <-    list(y = y[1:nspecies] , 
                      closure = visit$TP, site = visit$site_5km.n,
                      nclosure = nclosure, nsites = nsites, nspecies = nspecies,
                      SHORT = visit$SHORT, 
                      LONG = visit$LONG,
                      nobs = nrow(occup),
                      COVS=cov.array,
                      ncov=nmat)
  
# Specify model in BUGS language for vertical data format
mod = " model {
    ######################### PRIORS ##################################
    # ECOLOGICAL MODEL #
    # by species
    for (s in 1:nspecies){ 
      
      init.occ[s] ~ dunif(0,1) # initial occupancy
      
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
    
    tau.alpha.phi ~ dt(0,1,1)T(0,) # for SD truncated student t - approx half cauchy
    tau.gamma ~ dt(0,1,1)T(0,) 
    
    # covariate hyperpriors
    for(n in 1:ncov){
    beta.mean[n] ~ dunif(0,1) 
    mu.beta[n] <- logit(beta.mean[n]) # mean
    tau.beta[n] ~ dt(0,1,1)T(0,) # sd
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
    
    tau.p1 ~ dt(0,1,1)T(0,)
    tau.p2 ~ dt(0,1,1)T(0,)
    tau.p3 ~ dt(0,1,1)T(0,)
    
    tau.p ~ dt(0,1,1)T(0,) # for between year variance#
    
    ##################### END OF PRIORS #########################
    
    
    # Ecological model: Define state conditional on parameters
    for (s in 1:nspecies){ 
      for (i in 1:nsites){ # 5km site
        z[s,i,1] ~ dbern(init.occ[s])
        for (t in 2:nclosure){
          # persistence at i for t is function of pesticides
          logit(phi[s,i,t]) <- alpha.phi[s] + formula
          
          
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
  }  
  "
mod.sub <-gsub("formula",formula,mod)
cat(file = "model.txt",mod.sub)
cat(mod.sub)

# Initial values
zst <- zobs
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("phi","alpha.phi","mu.alpha.phi","tau.alpha.phi",
            "gamma","mu.gamma","tau.gamma",
            "alpha.p","tau.p",
            "beta","mu.beta","tau.beta","p",
            "init.occ")



source("rjags_func.R")


# Call JAGS and summarize posteriors
out <- jagsG(win.data, inits,n.adapt = na, params, "model.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb)

out$index=1

id <- paste0("_ID_",out$index)


if(daisy){
  s=1
  save(out,file=paste0(save.dir,"_C.",i,id,".rdata"))
}

save(out,file=rdataFile)

load("Jasmin_outputs/lad.B.all_C.2_run..rdata")
b=out
