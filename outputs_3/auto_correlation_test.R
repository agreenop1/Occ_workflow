################################################################################
################################################################################
# check autocorrelation
library(ncf)

source("combine_chain_4.1.R")

# parameters to keep
para <- c("init.occ","gamma")

# name of the start of the file
file="3_bee" 

# the number of iterations in each file
its <- c(500)

# the number of sites and number of time periods
time = 13
sites = 2278

# the number the iterations are thinned by
thin=5
red_it <- sum(its)/thin

# the filepath leading to the iterations
file.path = paste0("jas_out/",file ,"_C.CID_ID_UID.rdata")
file.path=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",file,"_C.CID_ID_UID.rdata")

# compile the chains from the model
out <- comb_daisy(parameters=para,
                  iter.index=10:10,chain.index=1:3,summary=F,file.path=file.path,by.it=NULL,
                  it.used=its,
                  iterations=NULL,verbose=T)
out <- simulation.l(out$samples)

z <- out$z

# gamma
gamma <- out$gamma

# initial occupancy
init <- out$init.occ
rm(out)

# phi
phi <- readRDS('phi.rds')

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
site_id <- read.csv("site_id.csv")
site_ll <- cbind(site=site_id$gr, BRCmap::gr2gps_latlon(site_id$gr))
write.csv(site_ll,'site_ll.csv')
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


# read out filepath
file.path=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/summary_outputs/summary.rds")
