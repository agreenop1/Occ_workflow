source("Occ_workflow/combine_chain_4.1.R")

# data formatted for jags
jas_data <- readRDS("Model_data/data_ladybird_all.499_1994.2010.rds")

file.path = paste0("Jasmin_outputs/", "lad.B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("beta","deviance","init.occ","gamma")
                  ,iter.index=1:28,chain.index=1:3,summary=T,file.path=file.path,by.it=500,iterations=4000,verbose=T)

# data in format from jags data prep
occup<- jas_data[[1]] #occ data
visit<- jas_data[[2]]# visit info
zobs<- jas_data[[3]]# init values
closure.period <- jas_data[[5]] #closure period

# check sample sizes in original data set
nsites <- length(unique(visit$site_5km)) # sites
nclosure <- length(unique(visit$TP)) # number of time points
nspecies <- ncol(occup[-1]) # num. species
nits <- out$mcmc.info$n.samples

out$model[[1]][[1]]$model

# covariate formula
covars <- jas_data[[4]]# covariates
covs <- list(temp=covars$temp,rqsum=covars$rq_sum_wide,semi=covars$semi) # named list
inters <- NULL # can supply a vector of which variables have interactions based on position in above list
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

psi <- array(dim=c(nits,nspecies,nsites,nclosure))

for (i in 1:nsites){
psi[,,i,1] <- out$sims.list$init.occ
}
phi<- array(dim=c(nits,nspecies,nsites,nclosure))
gamma <- out$sims.list$gamma
alpha.phi <- out$sims.list$gamma
beta <- out$sims.list$beta
COVS = cov.array

out <- NULL 

for (s in 1:nspecies){
 for (i in 1:nsites){
  for (t in 2:nclosure){
    phi[,s,i,t] <-   plogis (alpha.phi[,s] +
                          beta[,1,s]*COVS[i,t-1,1] +
                          beta[,2,s]*COVS[i,t-1,2] +
                          beta[,3,s]*COVS[i,t-1,3]) 
  
  psi[,s,i,t] <- psi[,s,i,t-1]*phi[,s,i,t] + (1-psi[,s,i,t-1])*gamma[,s]
 
    }
  } 
}

# percent or 0 beta
func.det <- array(dim=c(nclosure,nspecies,nsites,nits))
for(ni in 1:nits ){
  for (s in 1:nspecies){
    for (i in 1:nsites){
      for (t in 1:nclosure){
  
        func.det[t,s,i,ni] <- psi[ni,s,i,t]
  
      }
    }
  }
}
func.det[c(1,10),,1,1]
