################################################################################
# Summarise output of comb.chain function                                      #
################################################################################
summary_chains <- function(x,iterations="all",comb.chain=T,keep.samples='all'){
  
  n.burnin = x$n.burnin
  n.thin = x$thin
  n.chains = x$n.chains
  n.iter = x$n.iter
  DIC = T
  parameters.to.save = x$parameters.to.save
  codaOnly = F
  store.data = F
  samps <- x$samples
  
  
  
  if(iteration!="all"){
    
    samp.it = iteration/n.thin
    all.it = length(samps[[1]][,1,][[1]] )
      
    for (i in 1:n.chains){
      
    samps[[i]] <- samps[[i]][samp.it:all.it,,] 
    
    }
    cat("Iterations used...", length(samp.it:all.it))
  }
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