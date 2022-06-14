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

file="ag.lad" 
data="ag_ladybirds"

source("Occ_workflow/combine_chain_4.1.R")

file.path = paste0("Jasmin_outputs/",file ,".B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
                               "dtype1.p","dtype2.p","dtype3.p",
                               "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
                               "gamma","mu.gamma","tau.gamma",
                               "beta","mu.beta","tau.beta",
                               "init.occ"),
                  iter.index=1:21,chain.index=1:3,summary=T,file.path=file.path,by.it=1000,iterations=10000,verbose=T)

# output <- summary_chains(out,comb.chain = T,keep.samples = F)
options(max.print=10000)
# model 
print(out)
occ.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,]
samples = out$samples
samp.df <- ggmcmc::ggs(samples)

# co-variates in order used in model
covs <- c("Temperature_anom","Temp_mean","Semi.natural","Agri", "RQsum")

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2010.rds"))
species <- colnames(jas_data[[1]][-1])
nspecies <- length(species)
parameters=c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
             "dtype1.p","dtype2.p","dtype3.p",
             "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
             "gamma","mu.gamma","tau.gamma",
             "beta","mu.beta","tau.beta",
             "init.occ")



par.p <-sapply(parameters,function(x){
cat(x,"\n")
rhat <- out$Rhat[x][[1]]
par.plots <- list()

if(length(dim(rhat))>1){
 

  
  for (b in 1:length(covs)){
    
    cvr <- covs[[b]] 
    par.plots[[cvr]] <-list()
    plots <- list()
      
    for (i in 1:nspecies){
        
      par  <- paste0(x,"[",b,",",i,"]")
      tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(paste(species[i], "Rh",round(rhat[b,i],2)))+ theme(plot.title = element_text(size=9))
      ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)
      
      plots[[i]] <- ggpubr::ggarrange(tr,ds)
    }
    par.plots[[cvr]] <- ggpubr::ggarrange(plotlist = plots,nrow = 4)
  }
par.plots
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
    
    tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(paste("Rh =",r))+ theme(plot.title = element_text(size=9))
    ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
      geom_vline(xintercept = 0, linetype="dashed", 
                 color = "black", size=0.5)
    ct <-  ggs_caterpillar(samp.df[samp.df$Parameter==par,])+ 
      geom_vline(xintercept = 0, linetype="dashed", 
                 color = "black", size=0.5)
    plots[[i]] <- ggpubr::ggarrange(tr,ds,ct)
  }
  par.plots[[x]] <- if(par.l>1){plots}else{
    ggpubr::ggarrange(plotlist = plots,nrow = 4)}
  par.plots
  
} else if (length(rhat)==nspecies){
  

  par.plots[[x]] <-list()
  
  plots <- list()
  for (i in 1:nspecies){
    
    par  <- paste0(x,"[",i,"]")
    tr <- ggs_traceplot(samp.df[samp.df$Parameter==par,])+ggtitle(paste(species[i], "Rh",round(rhat[i],2)))+ theme(plot.title = element_text(size=9))
    ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
      geom_vline(xintercept = 0, linetype="dashed", 
                 color = "black", size=0.5)
    
    plots[[i]] <- ggpubr::ggarrange(tr,ds)
  }
  par.plots[[x]] <- ggpubr::ggarrange(plotlist = plots,nrow = 4)
} 
par.plots
},USE.NAMES=T)


#################################################################
############# Mean coefficient for each species #################
#################################################################


plot_effects <- function(x,species,cov){
  
  pr <-cov
  
  
  
  means <- out$mean$beta[pr,]
  lb <- out$q2.5$beta[pr,]
  ub <- out$q97.5$beta[pr,]
  
  pdat <- data.frame(matrix (nrow=length(species),ncol=4))
  colnames(pdat) <- c("Species","Mean","Lb","Ub")
  pdat[,1] <- species
  pdat[,2] <- means
  pdat[,3] <- lb
  pdat[,4] <- ub
  
  
  ggplot (pdat,aes(x=fct_rev((reorder(Species,Mean))),y=Mean))+ 
    geom_point(stat="identity",aes(x=fct_rev((reorder(pdat$Species,pdat$Mean)))))+
    geom_errorbar(aes(ymin=Lb,ymax=Ub)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),text=element_text(size=9) ,
          axis.text.x = element_text(angle = 90),axis.line = element_line(colour = "black"))+
    ylab(x)+xlab("Species")+ 
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "black", size=0.5)+coord_flip()
}

wid = 8
hei = 4
  
ta <- plot_effects("Temperature anomaly (deviation from mean of site)",species=species,cov=1)
ggsave(paste0("spec_plot/temp_a_",file,".png"),ta,width = wid,height = hei)

mt <- plot_effects("Mean temperature at a site",species=species,cov=2)
ggsave(paste0("spec_plot/mean_temp_",file,".png"),plot=mt,width = wid,height = hei)

sn <- plot_effects("Semi natural",species=species,cov=3)
ggsave(paste0("spec_plot/semi_nat_",file,".png"),plot=sn,width = wid,height = hei)

ag <- plot_effects("Agriculture",species=species,cov=4)
ggsave(paste0("spec_plot/agri_",file,".png"),plot=ag,width = wid,height = hei)

rq <- plot_effects("Sum risk quotient",species=species,cov=5)
ggsave(paste0("spec_plot/RQ_",file,".png"),plot=rq,width = wid,height = hei)



ggsave(paste0("spec_plot/mu_anom_",file,".png"),plot=par.p$mu.beta$mu.beta[[1]],
       width = wid,height = hei)

ggsave(paste0("spec_plot/mu_mean.t_",file,".png"),plot=par.p$mu.beta$mu.beta[[2]],
       width = wid,height = hei)

ggsave(paste0("spec_plot/mu_semi.nat_",file,".png"),plot=par.p$mu.beta$mu.beta[[3]],
       width = wid,height = hei)

ggsave(paste0("spec_plot/mu_agri.nat_",file,".png"),plot=par.p$mu.beta$mu.beta[[4]],
       width = wid,height = hei)

ggsave(paste0("spec_plot/mu_rq_",file,".png"),plot=par.p$mu.beta$mu.beta[[5]],
       width = wid,height = hei)

ggsave(paste0("spec_plot/mu_rq_",file,".png"),plot=par.p$mu.beta$mu.beta[[5]],
       width = wid,height = hei)