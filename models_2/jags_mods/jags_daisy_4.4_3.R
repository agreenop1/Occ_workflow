# A function for running occupancy models in serial updates or 'daisy chains'
################################################################################
# rdataFile - An R object containing a jags object, or the path to such
# an object (.rdata file containing an object called "out").  

i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

library(rjags)
library(R2jags)
set.seed(i)

sp.group <- c("bee",
              "spi",
              "hov",
              "lad")

group = sp.group[[4]]
mod = 4

# output files
rdataFile=paste0("Model_outputs/",mod,"_",group,"_all_C.",i,"_run",".rdata")
save.dir=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",mod,"_",group)



load(rdataFile) # file name

# run it = iterations each run; adapt = adaptive iterations; total it = max number of iterations; thinning factor
# save it = how often should iterations be saved - use to control memory + wall time (note model
# is not recompiled here)

by.it = 4000; adapt=100 ; total.it = 20000; n.thin = 5; save.it=100

# recompile the object
R2jags:::recompile.rjags(out,n.iter = adapt)
sufficient.adapt <- adapt(object=out$model,n.iter=1,end.adaptation=TRUE)
print(sufficient.adapt)
# Number run 
n.run = out$n.iter
out$n.iter

# Number left to run
n.left = total.it-n.run
n.left

# Number of runs left
n.left/by.it

# module needed
load.module("dic")

# iterations saved 
save.loop <- by.it/save.it # saves memory

for (s in 1:save.loop){
  
  # run index
  index <- out$index+1  
  
  if(out$n.iter<20000){
  
  # parameters
  params <- c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
              "dtype1.p","dtype2.p","dtype3.p",
              "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
              "gamma","mu.gamma","tau.gamma",
              "beta","mu.beta","tau.beta",
              "init.occ","deviance")
  }else {
    # parameters
    params <- c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
                "dtype1.p","dtype2.p","dtype3.p",
                "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
                "gamma","mu.gamma","tau.gamma",
                "beta","mu.beta","tau.beta",
                "init.occ","deviance","phi")  
  }
  
  # update the JAGS model by however many iterations
  samples <- coda.samples(out$model, variable.names=params, n.iter=save.it, thin = n.thin)
  
  #fit$n.burnin <- object$n.iter
  
  out <- list(model=out$model, n.thin=n.thin, parameters.to.save=params,
              model.file = out$model.file, n.iter=save.it+out$n.iter, DIC = out$DIC,samples=samples)
  
  
 
  
  class(out) <- "rjags"
  
  # index for saving to keep track of iterations
  out$index <- index
  out$adapt <- adapt
  out$sufficient.adapt <- sufficient.adapt
  
  id<- paste0("_ID_",out$index)
  
  # remove iterations count from file name
  out$iteration_id <- paste0(out$n.iter-by.it,"-",out$n.iter)
  
  
  save(out,file=paste0(save.dir,"_C.",i,id,".rdata"))
}

save(out,file=rdataFile)

cat("Number of iterations run =",out$n.iter, "Number of iterations left =",total.it- out$n.iter,
    " Number of runs left = ",(total.it- out$n.iter)/by.it)




