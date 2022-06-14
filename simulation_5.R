#########################################################################
# 5. Code for simulation of dynamics bayesian occupancy models          #
#########################################################################
library(BRCmap)
library(dplyr)
library(ggplot2)
library(jagsUI)
library(forcats)
library(dplyr)
library(ggmcmc)

# iter.index = UID in file path needs to be related to unique index containing different chunks of the chain (iteration numbers)
# chain.index = if comb.chain = T chains need to be combined. CID used to identify the different chain in the file path
# file.path = e.g. paste0("Model_outputs/Mod_sumRQ_500_C.","CID","_ID_","UID",".rdata")
# by.it = number of increases in iterations in each file
# iterations = iterations to be used
# verbose = T/F output to console

################################################################################
source("Occ_workflow/combine_chain_4.1.R")

file.path = paste0("Jasmin_outputs/", "lad.B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("beta","deviance","gamma","init.occ","alpha.phi"),
                  iter.index=1:13,chain.index=1:3,summary=T,file.path=file.path,by.it=1000,iterations=1000,verbose=T)

#  occ data
jas_data <- readRDS("Model_data/data_ladybirds_all.499_1994.2010.rds")

# output <- summary_chains(out,comb.chain = T,keep.samples = F)

# model 
print(out)

#############################################################
# Function to simulate occupancy under different covariates #
#############################################################

sim.occ <-  function(alpha.phi,psi,z,gamma,beta,n.beta,covars,nclosure,
                     samp.it=100,nspecies,nsites){
  
  
  COVS <- covars 
  phi <- array(dim=c(samp.it,nspecies,nsites,nclosure))
  
  for (s in 1:nspecies){
    for (i in 1:nsites){
      for (t in 2:nclosure){
        
        # intercept
        logit.phi <-  (alpha.phi[,s][1:samp.it]) 
        
        # covariate effects
        for (b in 1:n.beta){
          
          logit.phi <-  logit.phi + (beta[,b,s][1:samp.it])*COVS[i,t-1,b] 
        }
        
        # persistence
        phi[,s,i,t] <-   plogis (logit.phi) 
        
        # population occupancy
        psi[,s,i,t] <- psi[,s,i,t-1]*phi[,s,i,t] + (1-psi[,s,i,t-1])*(gamma[,s][1:samp.it]) 
        
        # realised occupancy 1 or 0
        prb  <- z[,s,i,t-1]*phi[,s,i,t] + (1-z[,s,i,t-1])*(gamma[,s][1:samp.it])
        
        z[,s,i,t] <- sapply(prb ,function(x){rbinom(1,1,x)})
      }
    } 
  }
  list(phi=phi,z=z,psi=psi)
}

#######################################################################################
# data in format from jags data prep
occup<- jas_data[[1]] #occ data
visit<- jas_data[[2]]# visit info
closure.period <- jas_data[[5]] #closure period
site_id <- jas_data[[6]]

# check sample sizes in original data set
nsites <- length(unique(visit$site_5km)) # sites
nclosure <- length(unique(visit$TP)) # number of time points
nspecies <- ncol(occup[-1]) # num. species
species <- colnames(occup[-1])
nits <- out$mcmc.info$n.samples

# covariate formula
covars <- jas_data[[4]]# covariates
covs <- list(covars$temp_anom,covars$mean_temp, covars$RQsum,covars$semi)
###################################################################
# covarite processing

# create array for covs
row = nrow(covars[[1]]); col = ncol(covs[[1]]); nmat = length(covs)
cov.array <- array(dim=c(row,col,nmat))

for(n in 1:nmat){
  cov.array[,,n]  <- as.matrix(covs[[n]])
}


###############################################################
# simulation params
sims <- out$sims.list # sims list
samp.it=60 # samples to use

# starting arrays
psi <- array(dim=c(samp.it,nspecies,nsites,nclosure))
z <- array(dim=c(samp.it,nspecies,nsites,nclosure))


# starting values
for (i in 1:nsites){
  psi[,,i,1] <- sims$init.occ[1:samp.it,]
  z[,,i,1] <- sapply(psi[,,i,1],function(x){rbinom(1,1,x)})
  
}

# simulate occupancy
base <- sim.occ(alpha.phi=sims$alpha.phi,psi=psi,z=z,gamma=sims$gamma,
                beta=sims$beta,n.beta=4,
                covars=cov.array,nclosure=10,
                samp.it=99,nspecies=nspecies,nsites=nsites)

phi.1=base$phi

phi.1[,,1,2]==phi.x[,,1,2]
################################################################
covar_adapt <- function(COVS,cov.change,rep.closure,ref.closure){
rep.col <- function(x,n){matrix(rep(x,each=n),ncol=n,byrow=T)}
                                
covs.new <-lapply(cov.change, function(x){ 
               
                    cov.i=x[[1]]
                    type=x[[2]]
                    per=x[[3]]
                
                             if(type=="constant"){
                               
                              var  =  COVS[,ref.closure,cov.i]
                              var.t = var*per 
                              mx = rep.col(var.t,rep.closure)
                              mx
                               
                             }else if(type=="step"){
                               
                              var  =  COVS[,ref.closure,cov.i]
                              ct = 1-per
                              ct1 = 1-(ct/rep.closure)
                              mx  = matrix(ncol = rep.closure,nrow=length(var))
                              mx[,1]<- var*ct1
                              for (i in 2:rep.closure){  mx[,i] <- mx[,i-1]*ct1} 
                               mx
                               
                             }else if(type=="selection"){
                               
                               mx  =  COVS[,per,cov.i]
                               mx
                               
                            }  else {
                              
                              mx  =  COVS[,ref.closure,cov.i]
                              mx[mx>-9999] = 0
                              mx = rep.col(mx,rep.closure)
                              mx}
                             })
         
  new.cov.array <- array(dim=c(row,rep.closure,nmat))
  for(n in 1:nmat){
   new.cov.array[,,n]  <- as.matrix(covs.new[[n]])
  }
  new.cov.array
} 
  
  if(!is.null(coef.change)){ lapply(coef.change,function(x){beta.new[,x[1],]  <- x[2]})}
###########################################################################################  
COV1<- covar_adapt(rep.closure=9, 
            cov.change= list(list(1,"step",1.05),list(2,"mean",1),list(3,"step",1.05),list(4,"constant",1)),
            ref.closure = 9, # as column reference 
            COVS = cov.array)
COV2<- covar_adapt(rep.closure=9, 
                   cov.change= list(list(1,"mean",1.10),list(2,"mean",1),list(3,"step",1.10),list(4,"constant",1)),
                   ref.closure = 9, # as column reference 
                   COVS = cov.array)

for (i in 1:nsites){
  psi[,,i,1] <- base$psi[,,i,10]
  z[,,i,1] <- sapply(psi[,,i,1],function(x){rbinom(1,1,x)})
  
}
occ1=sim.occ(alpha.phi=sims$alpha.phi,psi=psi,z=z,gamma=sims$gamma,
        beta=sims$beta,n.beta=4,
        covars=COV1,nclosure=10,
        samp.it=300,nspecies=nspecies,nsites=nsites)


occ2=sim.occ(alpha.phi=sims$alpha.phi,psi=psi,z=z,gamma=sims$gamma,
        beta=sims$beta,n.beta=4,
        covars=COV2,nclosure=10,
        samp.it=300,nspecies=nspecies,nsites=nsites)

#########################################################################
samples <-out$samples
nchains <- 3
samp.it <- length(samples[,"beta[1,1]",][[1]])
chain.par <- list()
COVS <- cov.array
n.beta <- 4  
nspecies = 2

  for(c in 1:nchains){
    
  phi <- array(dim=c(samp.it,nspecies,nsites,nclosure))
  psi <- array(dim=c(samp.it,nspecies,nsites,nclosure))
  z <- array(dim=c(samp.it,nspecies,nsites,nclosure))
  
  for (s in 1:nspecies){
    for (i in 1:nsites){
      psi[,s,i,1] <- samples[,paste0("init.occ[",s,"]"),][[c]]
      z[,s,i,1] <- sapply(psi[,s,i,1],function(x){rbinom(1,1,x)})
      
      for (t in 2:nclosure){
        
        # intercept
        logit.phi <-   samples[,paste0("alpha.phi[",s,"]"),][[c]]
        
        # covariate effects
        for (b in 1:n.beta){
          
          logit.phi <-  logit.phi + samples[,paste0("beta[",b,",",s,"]"),][[c]]*COVS[i,t-1,b] 
        }
        
        # persistence
        phi[,s,i,t] <-   plogis (logit.phi) 
        
        # population occupancy
        psi[,s,i,t] <- psi[,s,i,t-1]*phi[,s,i,t] + (1-psi[,s,i,t-1])*samples[,paste0("gamma[",s,"]"),][[c]] 
        
        # realised occupancy 1 or 0
        prb  <- z[,s,i,t-1]*phi[,s,i,t] + (1-z[,s,i,t-1])*samples[,paste0("gamma[",s,"]"),][[c]]
        
        z[,s,i,t] <- sapply(prb ,function(x){rbinom(1,1,x)})
      }
    } 
  }
 chain.par[[c]] <-list(z=z,phi=phi,psi=psi)
}
s1= samples[,"phi[1,1,2]",][[1]]
chain.par[[1]]$phi[,1,1,2]
dimnames(samples)


  logit.phi <-  samples[,paste0("alpha.phi[",s,"]"),][[c]] +
                samples[,paste0("beta[",1,",",s,"]"),][[c]]*COVS[i,t-1,1] +
                samples[,paste0("beta[",2,",",s,"]"),][[c]]*COVS[i,t-1,2] +
                samples[,paste0("beta[",3,",",s,"]"),][[c]]*COVS[i,t-1,3] +
                samples[,paste0("beta[",4,",",s,"]"),][[c]]*COVS[i,t-1,4] 

# persistence
s2= plogis(logit.phi) 
  i=1
  t=2
 s=1
 c=1
 
plot(density(s1)) 
plot(density(s2)) 
