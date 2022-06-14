########################################################################
# 4. Code for handling outputs of dynamics bayesian occupancy models  #
########################################################################
library(dplyr)
library(jagsUI)
library(forcats)
library(ncf)

source("combine_chain_4.1.R")
# iter.index = UID in file path needs to be related to unique index containing different chunks of the chain (iteration numbers)
# chain.index = if comb.chain = T chains need to be combined. CID used to identify the different chain in the file path
# file.path = e.g. paste0("Model_outputs/Mod_sumRQ_500_C.","CID","_ID_","UID",".rdata")
# by.it = number of increases in iterations in each file
# iterations = iterations to be used
# verbose = T/F output to console


# parameters to keep
para <- c("init.occ","gamma","beta","mu.beta")

# name of the start of the file
file="3_bee" 

# the filepath leading to the iterations

file.path=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",file,"_C.CID_ID_UID.rdata")


# the number of iterations in each file
its <- c(500)
iter.index <- 10:10
chains <- 3


out <- comb_daisy(parameters=para,
                  iter.index=iter.index,chain.index=1:chains,summary=T,file.path=file.path,by.it=NULL,
                  it.used=its,
                  iterations=NULL,verbose=T)

################################################################################
################################################################################
# the number of sites and number of time periods
time = 13
sites = 2278

# the number the iterations are thinned by
thin=5

#
red_it <- sum(its)/thin


# parameters to keep
para <- c("init.occ","gamma","z","phi")

# compile the chains from the model
out <- comb_daisy(parameters=para,
                  iter.index=iter.index,chain.index=1:chains,summary=F,file.path=file.path,by.it=NULL,
                  it.used=its,
                  iterations=NULL,verbose=T)
out <- simulation.l(out$samples)

# out latent occupancy
z <- out$z

# gamma
gamma <- out$gamma

# phi
phi <- out$phi

# initial occupancy
init <- out$init.occ
rm(out)


# setup psi  data
psi <- array(dim = dim(phi))

# initialize psi
for(i in 1:sites){psi[,,i,1] <- init}

# calculate population occupancy
for(t in 2:time){
  for(i in 1:sites){
    
    psi[,,i,t] <- psi[,,i,t-1]*phi[,,i,t] + (1-psi[,,i,t-1])*gamma 
  }
}


# calculate residuals
resid <- z - psi

# take residual mean
spacial.test <- apply(resid,c(1,3,4),mean)

# site coordinates
site_id <- read.csv("site_ll.csv")

# set up correlation matrix
max_d <- 680 # km 
increment_d <- 10 # number of distance bins (km)
nd <- max_d/increment_d
cmx <- array(dim=c(red_it,nd,time))
cmxr <- array(dim=c(red_it,nd,time))

for(i in 1:nd){
  for(j in 2:time){
    
    # test residuals
    
    spat.3 <- cbind(res=spacial.test[i,,j],site_id[3:4])
    
    # calculate MORAN I 
    corl <- correlog(x=spat.3$LONGITUDE,y=spat.3$LATITUDE,spat.3$res,latlon = T,increment = increment_d,resamp = 0)
    
    # shuffle iteration
    shuf <- spat.3$res[sample(1:length(spat.3$res),length(spat.3$res))]
    
    # calculate correlation on shuffle iteration
    corl_s <- correlog(x=spat.3$LONGITUDE,y=spat.3$LATITUDE,shuf,latlon = T,increment = increment_d,resamp = 0)
    
    
    cmx[i,1:68,j] <- corl$correlation
    cmxr[i,1:68,j] <- corl_s$correlation
  }
  
}

# assess temporal auto correlation
temporal.v <- aperm( spacial.test,c(3,2,1))[-1,,]
temporal.v <- aperm( temporal.v,c(2,3,1))
temp <-  apply(temporal.v,c(1,2),acf)

# shuffle values
temps <- apply(temporal.v,c(1,2),function(x){
  
  x = sample(x,replace = F)
  acf(x)
})


# read out filepath
file.path=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/summary_outputs/summary.rds")

# save output
saveRDS(list(
  
  out=out,
  spatial=list(observed=cmx,expected=cmxr),
  temporal=list(observed=temp,expected=temps)
)
        ,file.path)