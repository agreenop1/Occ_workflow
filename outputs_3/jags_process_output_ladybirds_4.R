########################################################################
# 4. Code for handling outputs of dynamics bayesian occupancy models  #
########################################################################
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

file="1_ladybirds" 
data="ladybirds"

source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")

file.path = paste0("jas_out/",file ,"_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
                               "dtype1.p","dtype2.p","dtype3.p",
                               "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
                               "gamma","mu.gamma","tau.gamma",
                               "beta","mu.beta","tau.beta",
                               "init.occ","deviance"),
                  iter.index=1:30,chain.index=1:3,summary=T,file.path=file.path,by.it=1000,iterations=15000,verbose=T)

# output <- summary_chains(out,comb.chain = T,keep.samples = F)
options(max.print=10000)
# model 
print(out)
samples = out$samples
occ.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,]
rhat. <- out$Rhat

samp.df <- ggmcmc::ggs(samples)

# co-variates in order used in model
if(file=="1_ladybirds"){
covs <- c("Mean Temperature","Temperature Anomaly"
          ,"Risk Quotient Anomaly","Risk Quotient Mean")
}else if(file=="2_ladybirds"){
  covs <- c("Mean Temperature","Temperature Anomaly",
            "Semi Natural Landcover"
            ,"Risk Quotient Anomaly","Risk Quotient Mean")
}else{
  covs <- c("Mean Temperature","Temperature Anomaly",
            "Semi Natural Landcover","Agricultural Landcover"
            ,"Risk Quotient Anomaly","Risk Quotient Mean"
            )
}

# number of covariates
n.cov = length(covs)  


# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))
species <- colnames(jas_data[[1]][-1])
nspecies <- length(species)
parameters=c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
             "dtype1.p","dtype2.p","dtype3.p",
             "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
             "gamma","mu.gamma","tau.gamma",
             "beta","mu.beta","tau.beta",
             "init.occ")


parameters=c("mu.beta","tau.beta","beta")

# par.p plots all parameters above - mainly diagnostics -  should detect each type of parameter okay
# parameters need to be a vector of names
par.p <-sapply(parameters,function(x,simple=T){
  cat(x,"\n")
  rhat <- out$Rhat[x][[1]]
  par.plots <- list()
  
  # beta coefs for species level parameters
  if(length(dim(rhat))>1){
    
    
    
    for (b in 1:length(covs)){
      
      cvr <- covs[[b]] 
      par.plots[[cvr]] <-list()
      plots <- list()
      
      for (i in 1:nspecies){
        
        # plots
        par  <- paste0(x,"[",b,",",i,"]")
        
        # traceplot
        tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(paste(species[i], "Rh",round(rhat[b,i],2)))+ theme(plot.title = element_text(size=9))
        
        # density plot
        ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
          geom_vline(xintercept = 0, linetype="dashed", 
                     color = "black", size=0.5)
        
        plots[[i]] <- ggpubr::ggarrange(tr,ds)
      }
      
      # arranged plots
      par.plots[[cvr]] <- ggpubr::ggarrange(plotlist = plots,nrow = 4)
    }
    par.plots
    
    # hyperparameters or those indexed higher than species
  }else if (length(rhat)<nspecies){
    par.l  <- length(rhat)
    par.plots[[x]] <-list()
    plots <- list()
    for (i in 1:par.l){
      
      if(par.l>1){ par  <- paste0(x,"[",i,"]")
      r <- round(rhat[i],3)
      
      }else{
        par <-x
        r   <- round(rhat,3)}
      
      # plots
      tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(covs[[i]],subtitle=paste("Rh =",r))+
        theme(plot.title = element_text(size=9),plot.subtitle = element_text(size=8))
      ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)
      ct <-  ggs_caterpillar(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)
      
      # arrange plots
      if(simple==F){
      plots[[i]] <- ggpubr::ggarrange(tr,ds,ct)
      }else{
        plots[[i]] <- ggpubr::ggarrange(tr,ct)
      }
    }
    par.plots[[x]] <- if(par.l>1){plots}else{
      ggpubr::ggarrange(plotlist = plots,nrow = 4)}
    par.plots
    
    # species level parameters
  } else if (length(rhat)==nspecies){
    
    rhat
    par.plots[[x]] <-list()
    
    plots <- list()
    for (i in 1:nspecies){
      
      # plots
      par  <- paste0(x,"[",i,"]")
      tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(paste(species[i], "Rh",round(rhat[i],2)))+ theme(plot.title = element_text(size=9))
      ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)
      
      # arrange plots
      plots[[i]] <- ggpubr::ggarrange(tr,ds)
    }
    par.plots[[x]] <- ggpubr::ggarrange(plotlist = plots,nrow = 4)
  } 
  par.plots
},USE.NAMES=T)

# save plots
wid = 8
hei = 8

m.beta

# plots of chains, density and parameter estimate 
m.beta <- ggpubr::ggarrange(plotlist = par.p$mu.beta$mu.beta,nrow=4)
ggsave(paste0("Main Effects",'_',file,".png"),m.beta,width = wid,height = hei)
ggsave(paste0("Main Effects.2",'_',file,".png"),m.beta[[2]],width = wid,height = hei)
m.beta
#################################################################
############# Mean coefficient for each species #################
#################################################################
# plot each species level coef

plot_effects <- function(x,species,cov){
  
  pr <-cov
  
  
  
  means <- out$mean$beta[pr,]
  lb <- out$q2.5$beta[pr,]
  ub <- out$q97.5$beta[pr,]
  rhats <- rhat.$beta[pr,]
  rhats <- ifelse(rhats>1.05,'not converged','converged')
  
  pdat <- data.frame(matrix (nrow=length(species),ncol=5))
  colnames(pdat) <- c("Species","Mean","Lb","Ub",'Convergence')
  pdat[,1] <- species
  pdat[,2] <- means
  pdat[,3] <- lb
  pdat[,4] <- ub
  pdat[,5] <- rhats
  
  
  ggplot (pdat,aes(x=fct_rev((reorder(Species,Mean))),y=Mean,color=Convergence))+ 
    geom_point(stat="identity",aes(x=fct_rev((reorder(pdat$Species,pdat$Mean)))))+
    geom_errorbar(aes(ymin=Lb,ymax=Ub)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),text=element_text(size=9) ,
          axis.text.x = element_text(angle = 90),axis.line = element_line(colour = "black"))+
    ylab(x)+xlab("Species")+ 
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "black", size=0.5)+coord_flip()
}




ind_plots <- list()

for(i in 1:n.cov){
ind_plots[[i]] <-  plot_effects(covs[[i]],species=species,cov=i)
ggsave(paste0(covs[[i]],'_',file,".png"),ind_plots[[i]],width = wid,height = hei)
}

ind_plots

## comparions plots for parameters in different models ##############################
x=out$sims.list; y = out$sims.list; parameter = "mu.beta"; samples = 500

npar <- ncol(as.data.frame(x["mu.beta"]))

if(!is.null(y)){
  sapply(1:5,function(x){
    
  })
}






