
# iter.index = UID in file path needs to be related to unique index containing different chunks of the chain (iteration numbers)
# chain.index = if comb.chain = T chains need to be combined. CID used to identify the different chain in the file path
# file.path = e.g. paste0("Model_outputs/Mod_sumRQ_500_C.","CID","_ID_","UID",".rdata")
# by.it = number of increases in iterations in each file
# thin = thin rate
# iterations = iterations to be used
# comb.chain = T/F combining chains run in parallel
# verbose = T/F output to console

###############################################################################
comb_daisy <- function(parameters=NULL,iter.index,chain.index,file.path,by.it,iterations="all",verbose){
  
  library("rjags")
  
  il <- length(iter.index) 
  
  
  
  # get total iterations from model
  
  # combining multiple chains that have been run in parallel
    c1.end <- gsub("UID",iter.index[il],file.path) # chain 1 last iteration
    c1.end <- gsub("CID",chain.index[1],c1.end)
    load(c1.end)
    t.iter <- out$n.iter # get total iterations 
    if(verbose){ cat("Total iterations =",t.iter,"...")}
    
    
  # do we use all or subset of iterations supplied in file index
  if(!iterations=="all"){
    
    seg=t.iter/by.it # how many objects are there
    keep.its =(iterations/by.it)-1 # how many should we keep
    if(keep.its==0){it.index=max(seg)}else{
    it.index=iter.index[(seg-keep.its):max(seg)]}
    if(verbose){cat("Using",iterations,"iterations...","\n")}
    
    
  }else{
    
    it.index=iter.index  
    if(verbose){cat("Using all iterations supplied in iter.index...","\n")}
    
  }
  
  # combine samples
  # variable select
 if(!is.null(parameters)){
      par <- paste0(parameters,collapse = "|")
      vars. <- varnames(out$samples)
      vars <- vars.[stringr::str_detect(vars.,par)]
      if(verbose){print(vars)}
  }
    
  # length of index accounting for iterations
  il <- length(it.index) 
  
  chains <- list()  # output for chains
    
  for (c in 1:length(chain.index)){
      
      mods  <- list()
      c.indx <- gsub("CID",chain.index[c],file.path) # set up chain file path
      if(verbose){cat("CHAIN",c,"\n")}
      
      for (i in 1:il){
        if(length(it.index==1)){file.name=gsub("UID",it.index,c.indx)}else{ # read in iteration
        file.name=gsub("UID",it.index[i],c.indx)} # read in iteration
        
        if (verbose){print(basename(file.name))}
        load(file.name)
        
        mods[[i]] <-  out # model object
        mods[[i]]$file.name <- file.name
        
         if(i==1){
          
          if(is.null(parameters)){ # store first samples
          samp <- mods[[i]]$samples }else{
          samp <- mods[[i]]$samples[,vars,drop=T] } 
          
          
          s.tart <- start(samp)
          
         }else{
          
          if(is.null(parameters)){ # add new samples else 
              samp <- jagsUI:::bind.mcmc(samp,mods[[i]]$samples,start=s.tart,n.new.iter=by.it)
              }else{ 
              samp <- jagsUI:::bind.mcmc(samp,mods[[i]]$samples[,vars,drop=T],start=s.tart,n.new.iter=by.it)}
              
           }
       }
      
     chains$mods[[c]] <- mods
     chains$samples[[c]] <- samp
     cat("Iterations used after thin", length(time(samp)),"\n")
    }
    
    if(is.null(parameters)){
    parameters.to.save <- mods[[1]]$parameters.to.save} else{
    parameters.to.save <- c(vars)
    }
    n.adapt=mods[[length(it.index)]]$adapt
    sufficient.adapt=mods[[length(it.index)]]$sufficient.adapt
  
    out<-list(models=chains$mods,samples=chains$samples,n.burnin=t.iter-iterations,n.iter=iterations,thin=thin(chains$samples[[1]]),
              n.chains=length(chain.index),parameters.to.save=parameters.to.save,n.adapt=n.adapt, sufficient.adapt= sufficient.adapt)
}    
######################################################
summary_chains <- function(x,comb.chain,keep.samples='all'){
  
  n.burnin = x$n.burnin
  n.thin = x$thin
  n.chains = x$n.chains
  n.iter = x$n.iter
  DIC = T
  parameters.to.save = x$parameters.to.save
  codaOnly = F
  store.data = F
  samps <- x$samples
  x$samples <- NULL
  mods <- x$models
  fin.mod = x$models[[1]][[length(x$models[[1]])]]
  s.adapt=fin.mod$sufficient.adapt
  t.adapt=100
  verbose=T
  
  #Create empty lists
  out <- samples <- model <- list()
  total.adapt <- sufficient.adapt <- vector(length=n.chains)
  
  if(comb.chain){  
    #Save samples and model objects from each cluster
    for (i in 1:n.chains){
      samples[[i]] <- samps[[i]][[1]]
      samps[[i]][[1]] <- NULL
      model[[i]] <- mods[[i]]
      total.adapt[i] <- t.adapt
      sufficient.adapt[i] <-  s.adapt
    }
    out$samples <- as.mcmc.list(samples)  
  }else{
    out$samples <-samps
    samps[[i]][[1]] <- NULL
    model=mods
    total.adapt <- t.adapt
    sufficient.adapt<-  s.adapt
  }
  # Remove columns with all NA
  try({
    all_na <- apply(as.matrix(out$samples),2, function(x) all(is.na(x)))
    out$samples <- out$samples[,!all_na]
  })
  out$model <- model
  out$total.adapt <- total.adapt
  out$sufficient.adapt <- sufficient.adapt
  names(out$model) <- sapply(1:length(out$model),function(i){paste('cluster',i,sep="")})
  
  
  samples <- out$samples
  m <- out$model
  total.adapt <- out$total.adapt
  sufficient.adapt <- out$sufficient.adapt
  
  
  #Get more info about MCMC run
  time <-  NULL
  
  
  #Combine mcmc info into list
  n.samples <- dim(samples[[1]])[1] * n.chains
  end.values <- samples[(n.samples/n.chains),]
  mcmc.info <- list(n.chains,n.adapt=total.adapt,sufficient.adapt,n.iter,n.burnin,n.thin,n.samples,end.values,time)
  names(mcmc.info) <- c('n.chains','n.adapt','sufficient.adapt','n.iter','n.burnin','n.thin','n.samples','end.values','elapsed.mins')
  
  #Reorganize JAGS output to match input parameter order
  if(dim(samples[[1]])[2]>1){
    samples <- jagsUI::: order.params(samples,parameters.to.save,DIC,verbose=verbose)
  }
  
  #Convert rjags output to jagsUI form 
  output <- jagsUI:::process.output(samples,DIC=DIC,codaOnly,verbose=verbose)
  if(is.null(output)){
    output <- list()
    samples <- order.params(samples,parameters.to.save,DIC,verbose=verbose)
    output$samples <- samples
    output$model <- m
    output$n.cores <- n.cores
    class(output) <- 'jagsUIbasic'
    return(output)
  }
  
  #Summary
  output$summary <- jagsUI:::summary.matrix(output,samples,n.chains,codaOnly)
  
  if(!is.null(keep.samples)){
        if(keep.samples=="all"){output$samples <- samples
     }else{output$samples <- samples}
    }
  output$modfile <- fin.mod$model.file
  
  #If user wants to save input data/inits
  if(store.data){
    output$inits <- inits
    output$data <- data
  } 
  output$model <- mods
  output$parameters <- parameters.to.save
  output$mcmc.info <- mcmc.info
  output$run.date <- NULL
  output$random.seed <- NULL
  output$parallel <- TRUE
  output$bugs.format <- F
  output$calc.DIC <- DIC
  mcmc.info$n.cores <- n.chains
  #Classify final output object
  class(output) <- 'jagsUI'
  
  return(output)
  
}



  
