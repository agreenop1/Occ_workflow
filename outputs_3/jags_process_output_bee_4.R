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


source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")

# read in parameters
para <- c("phi")

# start of file name
file="3_bee" 

# read in data with cover it information
data="bees"

# the number of sites and number of time periods
time = 13
sites = 2278

# number of iterations in each file
its <- c(500)

# file path where the files are located
file.path = paste0("jas_out/",file ,"_C.CID_ID_UID.rdata")

# iterations used
iterations. <- 200

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))

# compile the chains from the model
#out <- comb_daisy(parameters=para,
                  iter.index=20:20,chain.index=1:2,summary=F,file.path=file.path,by.it=500,
                  it.used=its,
                  iterations=1000,verbose=T)
out <- simulation.l(out$samples,5,2,out$parameters.to.save)

out <- readRDS("jas_out/summary_p.rds")[[1]]

# output <- summary_chains(out,comb.chain = T,keep.samples = F)
options(max.print=10000)
# model 
print(out)
samples = out$samples
eocc.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,]
rhat. <- out$Rhat

samp.df <- ggmcmc::ggs(samples)

# co-variates in order used in model
if(file=="1_bee"){
covs <- c("Mean Temperature","Temperature Anomaly"
          ,"Risk Quotient Anomaly","Risk Quotient Mean")
}else if(file=="2_bee"){
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


parameters=c("mu.beta")

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
      tr <- ggs_density(samp.df[samp.df$Parameter==par,])+ggtitle(covs[[i]],subtitle=paste("Rh =",r))+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)+
        theme(plot.title = element_text(size=9),plot.subtitle = element_text(size=8))+
        theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
      ds <- ggs_density(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)+
        theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
      ct <-  ggs_caterpillar(samp.df[samp.df$Parameter==par,])+ 
        geom_vline(xintercept = 0, linetype="dashed", 
                   color = "black", size=0.5)+
        theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
      
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

par.p$mu.beta$mu.beta

# plots of chains, density and parameter estimate 

m.beta <- ggpubr::ggarrange(plotlist = par.p$mu.beta.mu.beta,nrow=4)
ggsave(paste0("effect_plot/Main Effects.1",'_',file,".png"),m.beta[[1]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.2",'_',file,".png"),m.beta[[2]],width = wid,height = hei)
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
ggsave(paste0( "effect_plot/",covs[[i]],'_',file,".png"),ind_plots[[i]],width = wid,height = hei)
}

ind_plots

# check correlation between parameter estimates
cor_plot <- ggs_pairs(samp.df,family="mu.beta" ,lower = list(continuous = "density",alpha=0.2))
ggsave('check_plot/mu_cor_plots.png',cor_plot,width=12,height=9)

rm(out)
load("jas_out/3_bee_C.3_ID_9.rdata")
save(out,file="jas_out/3_bee_all_C.3_run.rdata")





################################################################################
################################################################################
# check autocorrelation
library(ncf)

source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")

para <- c("init.occ","gamma","beta","z")
file="3_bee" 
data="bees"
its <- c(500)
thin=5
red_it <- sum(its)/thin
file.path = paste0("jas_out/",file ,"_C.CID_ID_UID.rdata")


# compile the chains from the model
out <- comb_daisy(parameters=para,
                  iter.index=10:10,chain.index=1:1,summary=F,file.path=file.path,by.it=500,
                  it.used=its,
                  iterations=1000,verbose=T)
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

# the number of sites and number of time periods
time = 13
sites = 2278



# initialize psi
for(i in 1:sites){psi[,,i,1] <- init}

# calculate population occupancy
for(t in 2:time){
  for(i in 1:sites){
    
    psi[,,i,t] <- psi[,,i,t-1]*phi[,,i,t] + (1-psi[,,i,t-1])*gamma 
}
}


# calculate residuals
resid <- z - phi

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
cmx <- array(dim=c(red_it,nd))
cmxr <- array(dim=c(red_it,nd))

for(i in 1:nd){

# test residuals

spat.3 <- cbind(res=spacial.test[i,,-1],site_id[3:4])
x <- spacial.test[i,,-1]
# calculate MORAN I 
corl <- correlog(x=spat.3$LONGITUDE,y=spat.3$LATITUDE,x,latlon = T,increment = increment_d,resamp = 0)

# shuffle iteration
shuf <- x[sample(1:nrow(x),replace=F),]

# calculate correlation on shuffle iteration
corl_s <- correlog(x=spat.3$LONGITUDE,y=spat.3$LATITUDE,shuf,latlon = T,increment = increment_d,resamp = 0)


cmx[i,1:68] <- corl$correlation
cmxr[i,1:68] <- corl_s$correlation
  }
  


temporal.v <- aperm( spacial.test,c(3,2,1))[-1,,]
temporal.v <- aperm( temporal.v,c(2,3,1))

temporal.v <- temporal.v[1:4,,]
temp <-  apply(temporal.v,c(1,2),acf)

temps <- apply(temporal.v,c(1,2),function(x){
  
       x = sample(x,replace = F)
       acf(x)
       })
################################################################################
# plots 
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])


    
    sy <- x1[c("E","N",years_n[[i]])]
    
    
    colnames(sy)[1:3] <- c("E","N","Year")
    
    out[[i]] <- ggplot() +
      geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
      ylim(0, 700000)  + geom_tile(data = sy, 
                                   aes(x = E, y = N, fill =Year))+
      scale_fill_continuous(type = "viridis", name = years_n[[i]])

################################################################################
sim.occ <-  function(alpha.phi,psi,z,gamma,beta,n.beta,covars,nclosure,
                     samp.it=100,nspecies,nsites){
  
  
  COVS <- covars 
  phi <- array(dim=c(samp.it,nspecies,nsites,nclosure))
  
  
  
  
  for (s in 1:nspecies){ # number of species
    for (i in 1:nsites){ # number of sites
      for (t in 2:nclosure){ # number of closure periods
        
        # persistence intercept
        logit.phi <-  sample(alpha.phi[,s],samp.it) 
        
        # covariate effects
        for (b in 1:n.beta){
          
          logit.phi <-  logit.phi + sample(beta[,b,s],samp.it)*COVS[i,t-1,b] 
        }
        
        # persistence
        phi[,s,i,t] <-   plogis (logit.phi) 
        
        # population occupancy
        psi[,s,i,t] <- psi[,s,i,t-1]*phi[,s,i,t] + (1-psi[,s,i,t-1])*sample(gamma[,s],samp.it) 
        
        # realised occupancy 1 or 0
        prb  <- z[,s,i,t-1]*phi[,s,i,t] + (1-z[,s,i,t-1])*sample(gamma[,s],samp.it)
        
        z[,s,i,t] <- sapply(prb ,function(x){rbinom(1,1,x)})
      }
    } 
  }
  list(phi=phi,z=z,psi=psi)
}