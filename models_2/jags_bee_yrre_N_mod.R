############################################################################################
# 2. Multi-species occupancy model for risk quotient & environment with random species     #
# effect                                                                                   #
############################################################################################
rm(list = ls())
suppressPackageStartupMessages(library(rjags,warn.conflicts = F,quietly=T))
suppressPackageStartupMessages(library(R2jags,warn.conflicts = F,quietly=T))
# slurm input
i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

set.seed(i) # 1:3 depending on chain

group = 'bee_yrre_N'
mod = 3
cat(mod,group,"\n")

# data formatted for jags

  
jas_data <- readRDS("data_bees_all.499_1994.2016.rds")


# output files
rdataFile=paste0("Model_outputs/",mod,"_",group,"_all_C.",i,"_run",".rdata")
save.dir=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",mod,"_",group)

# CHECK MCMC PARAMETERS !! # 
# MCMC settings
ni <- 1000; nt <- 5 ; nb <- 0 ; nc <- 1; na <- 100; n.species="all"; daisy=T # can be all

# covariate formula
covars <- jas_data[[4]]# covariates

base <-list(list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,RQA=covars$RQsum_A,RQM=covars$RQsum_M),
            list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,RQA=covars$RQsum_A,RQM=covars$RQsum_M),
            list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M))

covs <- c(base[[mod]]) 
inters <- NULL # can supply a vector of which variables have interactions based on position in above list
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
row = nrow(covars[[1]]); col = ncol(covars[[1]]); nmat = length(covs)
cov.array <- array(dim=c(row,col,nmat))

for(n in 1:nmat){
  cov.array[,,n]  <- as.matrix(covs[[n]])
}

# formula for main effects
formula <- paste0("beta[",1:nmat,",s]*COVS[i,t-1,",1:nmat,"]",collapse ="+")

# if any interactions add to formula

if(!is.null(inters)){
  coms <- combn(inters,2)  
  
  for (n in 1:ncol(coms)){
    ints1 <- paste0("+beta[",n+nmat,",s]*COVS[i,t-1,",coms[1,n],"]*COVS[i,t-1,",coms[2,n],"]")
    formula <- paste0(formula,ints1) }
  
  nmat <- ncol(coms)+nmat
}

vars <- strsplit(formula,split="+",fixed=T)[[1]]

# output formula
cat("Ecological model covariate structure","\n")
cat(vars,sep="+ \n")
cat(formula,"\n")


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
                    ncov=nmat,
                    Year=visit$Year,
                    nYear=jas_data$nYear,
                    region=jas_data$region,
                    nregion=jas_data$nregion)

# Specify model in BUGS language for vertical data format
mod = " model {
    ######################### PRIORS ##################################
    # ECOLOGICAL MODEL #
    # by species
    for (s in 1:nspecies){ 
      
      init.occ[s] ~ dnorm(mu.alpha.psi,tau.alpha.psi) # initial occupancy
      
      logitgamma[s] ~ dnorm(mu.gamma,tau.gamma)  #  colonisation prior -heterogeneity between species
      logit(gamma[s]) <- logitgamma[s]
      
      alpha.phi[s] ~ dnorm(mu.alpha.phi,tau.alpha.phi) # intercept for species persistence model logit
      
      for (n in 1:ncov){  
       beta[n,s] ~ dnorm(mu.beta[n],tau.beta[n]) # covariate array prior
      }
      
    }
    
    # ecology model hyperpriors
    mu.alpha.phi ~ dnorm(0,0.01) # persistence intercept (expected value=0, so persistence=0.5)
    mu.alpha.psi ~ dnorm(0,0.01)
    mu.gamma ~ dnorm(0,0.01) # colonisation intercept
    
    tau.alpha.phi ~ dt(0,1,1)T(0,) # for SD truncated student t - approx half cauchy
    tau.alpha.psi ~ dt(0,1,1)T(0,)
    tau.gamma ~ dt(0,1,1)T(0,) 
    tau.region.psi ~ dt(0,1,1)T(0,) 
    
    # covariate hyperpriors
    for(n in 1:ncov){
    beta.mean[n] ~ dunif(0,1) 
    mu.beta[n] <- logit(beta.mean[n]) # mean
    tau.beta[n] ~ dt(0,1,1)T(0,) # sd
    }
    
      # random region effect - sum to zero
      for (t in 1:(nregion-1)) {
      region.psi[t] ~ dnorm(0,tau.region.psi)
    }
    
      region.psi[nregion] <- -sum(region.psi[1:(nregion-1)]) 
    
    # OBSERVATION MODEL
    # Observation model priors
    for(i in 1:nspecies){
      dtype1.p[i] ~ dnorm(mu.d1.p,tau.p1)
      dtype2.p[i] ~ dnorm(mu.d2.p,tau.p2)
      dtype3.p[i] ~ dnorm(mu.d3.p,tau.p3)
    }
    
    # random year effect - sum to zero
      for (t in 1:(nYear-1)) {
      alpha.p[t] ~ dnorm(0,tau.p)
    }
    
    alpha.p[nYear] <- -sum(alpha.p[1:(nYear-1)]) 
    
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
         
        logit(init.psi[s,i])   <- init.occ[s] + region.psi[region[i]]
        
        z[s,i,1] ~ dbern(init.psi[s,i])
        
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
        logit(p[s,i]) <- alpha.p[Year[i]] + dtype1.p[s] + dtype2.p[s]*SHORT[i] +
          dtype3.p[s]*LONG[i]
        
        y[i,s] ~ dbern(z[s,site[i],closure[i]]*p[s,i])
      }
    }
  }  
  "
mod.sub <-gsub("formula",formula,mod)
cat(file = paste0(group,nmat,'.txt'),mod.sub)
cat(mod.sub)

# Initial values
zst <- zobs
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
            "dtype1.p","dtype2.p","dtype3.p",
            "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
            "gamma","mu.gamma","tau.gamma",
            "beta","mu.beta","tau.beta",
            "init.occ")


# adapted rjags which separates adapt and  burn in
source("jags_rjags_func_2.1.R")


# Call JAGS and summarize posteriors
out <- jagsG(win.data, inits,n.adapt = na, params,paste0(group,nmat,'.txt'), n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, jags.seed = i)


out$index=1

# daisy ID
id <- paste0("_ID_",out$index)

# save daisy
if(daisy){
  save(out,file=paste0(save.dir,"_C.",i,id,".rdata"))
}

# save runtime
save(out,file=rdataFile)
