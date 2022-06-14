################################################################################
# Summarise output of comb.chain function                                      #
################################################################################
summary_chains <- function(x,iterations="all",comb.chain=T,keep.samples='all'){
  process.output. <- function(x,DIC,params.omit=NULL,verbose=TRUE) {
    
    out <- tryCatch({
      
      if(verbose){cat('Calculating statistics.......','\n')}  
      
      # Get parameter names
      params <- colnames(x[[1]])
      
      #Get number of chains
      m <- length(x)
      
      #Collapse mcmc.lists into matrix
      mat = do.call(rbind,x)
      
      #Get # of iterations / chain
      n <- dim(mat)[1] / m
      
      #Get parameter dimensions
      dim <-jagsUI::: get.dim(params)
      
      #Create new parameter name vectors to handle non-scalar params
      expand <- sapply(strsplit(params, "\\["), "[", 1)
      params.simple <- unique(sapply(strsplit(params, "\\["), "[", 1))
      
      #Functions for statistics
      qs <- function(x,y){as.numeric(quantile(x,y))}
      #Overlap 0 function
      ov <- function(x){findInterval(0,sort(c(qs(x,0.025),qs(x,0.975))))==1}
      #f function (proportion of posterior with same sign as mean)
      gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}}
      #n.eff function
      calcneff <- function(x,n,m){
        xp <- matrix(x,nrow=n,ncol=m)
        xdot <- apply(xp,2,mean)
        s2 <- apply(xp,2,var)
        W <- mean(s2)
        
        #Non-degenerate case
        if (is.na(W)){
          n.eff <- NA
        } else if ((W > 1.e-8) && (m > 1)) {
          B <- n*var(xdot)
          sig2hat <- ((n-1)*W + B)/n      
          n.eff <- round(m*n*min(sig2hat/B,1),0)
          #Degenerate case
        } else {
          n.eff <- 1
        }
        n.eff
      }
      
      #Gelman diag function
      gd <- function(i,hold){
        r <- try(gelman.diag(hold[,i], autoburnin=FALSE)$psrf[1], silent=TRUE)
        if(inherits(r, "try-error") || !is.finite(r)) {
          r <- NA
        }
        return(r)
      }
      
      #Make blank lists
      sims.list <- means <- rhat <- n.eff <- se <- as.list(rep(NA,length(params.simple)))
      q2.5 <- q25 <- q50 <- q75 <- q97.5 <- overlap0 <- f <- as.list(rep(NA,length(params.simple)))
      names(sims.list) <- names(means) <- names(rhat) <- names(n.eff) <- params.simple
      names(se) <- names(q2.5) <- names(q25) <- names(q50) <- names(q75) <- names(q97.5) <- params.simple
      names(overlap0) <- names(f) <- params.simple
      
      #This function modifies objects in global environment (output is discarded)
      #Calculates statistics for each parameter
      calc.stats <- function(i){
        
        #If parameter is not a scalar (e.g. vector/array)
        if(!is.na(dim[i][1])){
          
          #Get all samples
          sims.list[[i]] <<- mat[,expand==i,drop=FALSE]
          
          #if every iteration is NA, don't do anything else
          if(all(is.na(sims.list[[i]]))){return(NA)}
          
          #If more than 1 chain, calculate rhat 
          #Done separately for each element of non-scalar parameter to avoid errors
          if(m > 1 && (!i%in%params.omit)){
            hold <- x[,expand==i,drop=FALSE]
            nelements <- sum(expand==i)
            rhat.vals <- sapply(1:nelements,gd,hold=hold)
            names(rhat.vals) <- colnames(hold[[1]])
            rhat[[i]] <<- jagsUI:::populate(rhat.vals,dim[[i]])
          } else if (m == 1){
            hold <- x[,expand==i]
            rhat[[i]] <<- array(NA,dim=dim[[i]])
          }
          
          #Calculate other statistics
          ld <- length(dim(sims.list[[i]]))
          means[[i]] <<- jagsUI:::populate(colMeans(sims.list[[i]]),dim[[i]])
          if(!i%in%params.omit){
            se[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),sd),dim=dim[[i]])
            q2.5[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),qs,0.025),dim=dim[[i]])
            q25[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),qs,0.25),dim=dim[[i]])
            q50[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),qs,0.5),dim=dim[[i]])
            q75[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),qs,0.75),dim=dim[[i]])
            q97.5[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),qs,0.975),dim=dim[[i]])
            overlap0[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),ov),dim=dim[[i]])
            f[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),gf),dim=dim[[i]])
            n.eff[[i]] <<- jagsUI:::populate(apply(sims.list[[i]],c(2:ld),calcneff,n,m),dim=dim[[i]])   
          }
          
          sims.list[[i]] <<- jagsUI:::populate(sims.list[[i]],dim=dim[[i]],simslist=T,samples=dim(mat)[1])
          
          #If parameter is a scalar
        } else {
          
          if(m > 1 && (!i%in%params.omit)){rhat[[i]] <<- gelman.diag(x[,i],autoburnin=FALSE)$psrf[1]}
          
          sims.list[[i]] <<- mat[,i]
          
          if(all(is.na(sims.list[[i]]))){return(NA)}
          
          means[[i]] <<- mean(sims.list[[i]])
          if(!i%in%params.omit){
            se[[i]] <<- sd(sims.list[[i]])
            q2.5[[i]] <<- qs(sims.list[[i]],0.025)
            q25[[i]] <<- qs(sims.list[[i]],0.25)
            q50[[i]] <<- qs(sims.list[[i]],0.5)
            q75[[i]] <<- qs(sims.list[[i]],0.75)
            q97.5[[i]] <<- qs(sims.list[[i]],0.975)
            overlap0[[i]] <<- ov(sims.list[[i]])
            f[[i]] <<- gf(sims.list[[i]])
            n.eff[[i]] <<- calcneff(sims.list[[i]],n,m)}
        }
        
      }
      
      #Actually run function(nullout not used for anything)
      nullout <- sapply(params.simple,calc.stats)
      
      #Warn user if at least one Rhat value was NA
      rhat.sub <- unlist(rhat)[!is.na(unlist(means))]
      if(NA%in%rhat.sub&&verbose){
        options(warn=1)
        warning('At least one Rhat value could not be calculated.')
        options(warn=0,error=NULL)
      }
      
      #Do DIC/pD calculations if requested by user
      if(DIC & 'deviance' %in% params){
        dev <- matrix(data=mat[,'deviance'],ncol=m,nrow=n)   
        pd <- numeric(m)
        dic <- numeric(m)    
        for (i in 1:m){
          pd[i] <- var(dev[,i])/2
          dic[i] <- mean(dev[,i]) + pd[i]
        }    
        pd <- mean(pd)
        dic <- mean(dic)
        
        #Return this list if DIC/pD requested
        if(verbose){cat('\nDone.','\n')}
        return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
                    f=f,Rhat=rhat,n.eff=n.eff,pD=pd,DIC=dic))
      } else {
        #Otherwise return list without pD/DIC
        if(verbose){cat('\nDone.','\n')}
        return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
                    f=f,Rhat=rhat,n.eff=n.eff))
      }
      
    }, error = function(cond){
      message('Calculating statistics failed with the following error:')
      message(cond)
      message('\nOutput falling back to class jagsUIbasic\n')
      return(NULL)
    }
    )
    return(out)
  }

  require(jagsUI)
  n.burnin = x$n.burnin
  n.thin = x$thin
  n.chains = x$n.chains
  n.iter = x$n.iter
  DIC = T
  parameters.to.save = x$parameters.to.save
  codaOnly = F
  store.data = F
  samps <- x$samples
  
  
  
  if(iterations!="all"){
    
    samp.it = iterations/n.thin
    end.it = length(samps[[1]][,1,][[1]] )
    start.it = end.it-samp.it 
    for (i in 1:n.chains){
      
    samps[[i]]<- as.mcmc.list(mcmc(samps[[i]][[1]],start = start.it,end = end.it,thin=5))
    
    }
    cat("Iterations used after thin...", samp.it,"\n")
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
    samples <- jagsUI:::order.params(samples,parameters.to.save,DIC,verbose=verbose)
  }
  
  #Convert rjags output to jagsUI form 
  output <- process.output.(samples,DIC=DIC,codaOnly,verbose=verbose)
  if(is.null(output)){
    output <- list()
    samples <- jagsUI::: order.params(samples,parameters.to.save,DIC,verbose=verbose)
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
