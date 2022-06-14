# A function for running occupancy models in serial updates or 'daisy chains'
################################################################################
# rdataFile - An R object containing a jags object, or the path to such
# an object (.rdata file containing an object called "out").  

i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

library(nimble)
source('nimble_save_f.R')
set.seed(i)


group = 'hoverflies_pol'
mod = "quad"

# output files
rdataFile=paste0("Model_outputs/",mod,"_",group,"_all_C.",i,"_run",".rds")
save.dir=paste0("/gws/nopw/j04/ceh_generic/arrgre/sumRQenv/",mod,"_",group)



out <- readRDS(rdataFile) # file name

out <- list(parameters=params,modState=stateList,iterations=ni,modData=win.data,inits=inits,
            chain=i,id=id,nimbleModel=modx)

# recompile nimble
mod. <- nimbleModel(code = out$nimbleModel, name = "modx",  constants = out$modData,
                    inits = out$inits,returnDef = F)


# build mcmc - in this case all the default samplers are used
conmod <- configureMCMC(mod.,monitors=out$parameters,thin=out$thin)


# build mcmc - in this case all the default samplers are used
modmcmc <- buildMCMC(conmod)

# compile nimble model
cmod <- compileNimble(mod.)

# compile mcmc algorithm
cmcmc <- compileNimble(modmcmc,project = mod.)

## load the saved "state" file
stateList <- out$modState

modelState <- stateList$modelState
mcmcState <- stateList$mcmcState

## restore the saved "state" into the new model and new MCMC
setModelState(cmod, modelState)
setMCMCstate(conmod,cmcmc, mcmcState)

# run it = iterations each run; adapt = adaptive iterations; total it = max number of iterations; thinning factor
# save it = how often should iterations be saved - use to control memory + wall time (note model
# is not recompiled here)

by.it = 5000; n.thin = 5; save.it=500

# number of iterations run
out$iterations


# iterations saved 
save.loop <- by.it/save.it # saves memory

for (s in 1:save.loop){
  
  
  # run the algorithm 
  cmcmc$run(niter=save.it,reset=F,resetMV = T)
  
  # run index
  index <- out$index+1  
  
  # save the state of the model
  stateList <- list(modelState = getModelState(cmod),
                    mcmcState = getMCMCstate(conmod, cmcmc))  
  
  
  # daisy ID
  id <- paste0("_ID_",1)
  
  # remove samples
  samples=as.matrix(cmcmc$mvSamples)
  
  # index for saving to keep track of iterations

  
  id<- paste0("_ID_",index)
  
  
  # save model data
  out$modState <- stateList
  out$iterations <- out$iterations+save.it
  out$id <- id
  out$index <- index


  saveRDS(out,file=paste0(index,"out.rds"))
  saveRDS(samples,file=paste0(index,"samples.rds"))
  #saveRDS(out,file=paste0(save.dir,"_C.",i,id,".rds"))
  #saveRDS(samples,file=paste0(save.dir,"_C.",i,id,"_samples_",".rds"))
  #saveRDS(out,file=rdataFile)
}







