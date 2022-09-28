################################################################################
################ Process Bayesian Model Outputs ################################
################################################################################
# This script processes the outputs from the models and carries out a number of#
# posterior predictive checks and also plot some figures of the results. It is #
# worth noting that the size of some of the files can crash the memory so keep #
# an eye on this. The input the script takes is an output file from jags that  #
# that combines and summarizes the occupancy outputs.                          #


# needed packages
library(brms)
library(BRCmap)
library(ggpubr)
library(ggmcmc)
library(ggplot2)
library(dplyr)
library(forcats)
library(caret)

# these scripts contain some miscellaneous functions to help with summarizing the outputs
source("output_functions.R")
source("Occ_workflow_V2/outputs_3/combine_chain_4.2.R")
UK <-  readRDS("UK_map.rds")

# summary function
summary_f <- function(x){
  mean.x <- apply(x,1,mean)
  ql <- apply(x,1,quantile,0.025)
  qu <- apply(x,1,quantile,0.975)
  cbind(mean.x,ql,qu)}


## read in data ##
# output from jags  
occ_output <- readRDS(paste0("jas_out/3_bee_yrre_summary.rds"))#
out <- occ_output

# input file for jags
jas_data <- readRDS(paste0("Model_data/data_bees_all.499_1994.2016.rds"))



# remove id column
occ <- jas_data[[1]][-1] # observed occupancy
vis <- jas_data[[2]] # observed visit data 
occdat <- cbind(occ,vis) 


# the number of sites, time periods and species
time = 13
sites = nrow(jas_data[[4]]$temp_anom)
nspecies = ncol(occ)

# covariates in the order they were in the jags model (used in simulations)
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)


# number of iterations found in the files before thin
cat((nrow(out$sims.list$alpha.phi)*5)/3,"iterations used before thin per chain","\n")

################################################################################
########################### Diagnostic Plots ###################################
# In this section we run diagnostic plots on the chains to check convergence.  #

# species names 
species <- colnames(occ)


# names of covariates in the order they are in the model
covs <- c("Temperature spatial","Temperature temporal",
          "Semi Natural Landcover","Agricultural Landcover"
          ,"Risk Quotient temporal","Risk Quotient spatial"
)

# number of covariates
n.cov <- length(covs)

# options(max.print=10000);print(out)


# model samples
samples = out$samples

# summary of model samples
occ.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,] # check parameter convergence
rhat <- out$Rhat # all rhats from jags outputs

samp.df <- ggmcmc::ggs(samples) # all samples in a format that can be read by the plot function (see ggmcmc package)

### check convergence and chain diagnostics of required parameters ###
# parameters to check
parameters=c("alpha.phi", "mu.beta","init.occ" ) #, "region.psi"



## output plots ##
# x = parameter name  (parameters need to be a vector of names)
# samp.df = samples in ggmcmc format
# rhat = rhat from jags output 
# species_names=NULL
# simple = T/F determines where the density plot is included (not advised for species level parameters)

# ecological model parameters
par.p <- sapply(parameters,output_plots,samp.df = samp.df,rhat=rhat,species_names = species,simple = T,simplify = F)

ggarrange(plotlist=par.p$mu.beta,nrow=4,ncol=1) # mu beta
ggarrange(plotlist=sapply(par.p$beta,function(x) x[[6]],simplify = F),nrow=4,ncol=1) # beta
ggarrange(plotlist= par.p$alpha.phi ,nrow=4,ncol=1) # year intercept obs
ggarrange( plotlist=par.p$region.psi,ncol=1, nrow=4)  # region
ggarrange( plotlist=par.p$init.occ,ncol=1, nrow=4)  # initial occupancy

# observation model parameters
obs.p <- sapply(c("alpha.p","dtype1.p"),
                        output_plots,
                        samp.df =samp.df ,
                        rhat= rhat,species_names = species,simple =T,simplify = F)


ggarrange( plotlist=obs.p$dtype1.p,ncol=1, nrow=4) # dtype.1


## save plots ##
wid = 6
hei = 3

# plots of chains, density and parameter estimate 
m.beta <- ggpubr::ggarrange(plotlist = par.p$mu.beta,nrow=2,ncol=1)
ggsave(paste0("effect_plot/Main Effects.1",'_',tribe,'_',group,".png"),m.beta[[1]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.2",'_',tribe,'_',group,".png"),m.beta[[2]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.3",'_',tribe,'_',group,".png"),m.beta[[3]],width = wid,height = hei)
m.beta

## individual species plots ##
ind_plots <- list() # output list for species plots

# plot_effects is a function that plots the individual species parameter estimates for the covariates #
# covs = list of covariate names                                                                      #
# species = species names (in the order they are in the occupancy data)                               #
# out = jags output                                                                                   #
# rhat = rhat from jags output                                                                        #

# cycle through covariates and get individual species parameters
for(i in 1:n.cov){
  
  # output plot for species for each covariate
  ind_plots[[i]] <-  plot_effects(covs[[i]],species=species,cov=i,out = out,rhat=rhat)
  
  # save this output in a file called effect_plot
  ggsave(paste0( "effect_plot/",covs[[i]],'_',file,".png"),ind_plots[[i]],width = wid,height = hei)
}

# look at all individual species parameters
ind_plots

# check correlation between parameter estimates
cor_plot <- ggs_pairs(samp.df,family=c("mu.beta") ,lower = list(continuous = "density",alpha=0.2))
ggsave('check_plot/mu_cor_plots.png',cor_plot,width=12,height=9)


#################################################################
############# Posterior Predictive Checks #######################
# In this section we run some simulations to look at posterior  #
# predictive checks.                                            #

covs <- c(base) # covariate list

# the function we run here predicts occupancy based on the model parameters
# keep an eye on the number of samples (nrep) due to memory issues

#  number of samples
nrep <- 200 #dim(out$sims.list$beta)[1]

# parameters from the model
beta <- out$sims.list$beta[1:nrep,,] # beta coefficient (persistence)
gamma <-  out$sims.list$gamma[1:nrep,] # colonization
init <-  out$sims.list$init.occ[1:nrep,] # initial occupancy
alpha <- out$sims.list$alpha.phi[1:nrep,] # intercept (persistence) 
region.psi <- out$sims.list$region.psi[1:nrep,] # region (initial occupancy)
region.cov <- jas_data$nregion # number of regions

## this function runs the population simulations ##
# estimates of population occupancy for species


psi <-  sim.pop(psi.start = init,
                gamma = gamma,
                beta = beta,
                alpha = alpha,
                psi =  array(dim=c(nrep,nspecies,sites,time)),
                cov.array = f_array(covs),
                species = T,
                region = region.psi,
                region.cov = jas_data$region
               )



# predict occupancy state - write out occupancies to save memory (this is all a bit slow)
dir.create("binary_occupancy")
for(i in 1:nrep){
  
a.x  <-  apply(psi[i,,,],c(2,3),function(x){rbinom(length(x),1,x)})
saveRDS(a.x,file=paste0("binary_occupancy/",data,"_rep_",i,"_",mod,".rds"))

}

# observed occupancy
z <- jas_data[[3]]


# observation model

alpha.p <- out$sims.list$alpha.p
d1 <- out$sims.list$dtype1.p
d2 <- out$sims.list$dtype2.p
d3 <- out$sims.list$dtype3.p


# observation model parameters
observed <- jas_data[[2]] # visit data
SHORT <- observed$SHORT # list length for each visit
LONG <- observed$LONG
nobs <- length(SHORT) # number of observations
site <- observed$site_5km.n 
closure <- observed$TP
Year <- observed$TP
y <- array(dim=c(nobs,nspecies),0) # predicted observations
p <- array(dim=c(nobs,nspecies)) # probability occupancy


# aggregate all observations seen at a site
vis_rep <- list()
sum_rep <- list()
spe_rep <- list()

v_x_o <- vector()
v_x_e <- vector()
s_x_o <- vector()
s_x_e <- vector()
e <- 0.0001




# create directory for outputs
#dir.create("predicted_observed_occupancy")

# observation model
for(i in 1:nrep){
  
  zrep <- readRDS(paste0("binary_occupancy/",data,"_rep_",i,"_",mod,".rds"))
  
  for(o in 1:nobs){
    
    
    # observation probability  
    p[o,] <- inv_logit_scaled( alpha.p[i,Year[o]] + 
                                 d1[i,] + 
                                 d2[i,]*SHORT[o] +
                                 d3[i,]*LONG[o])
    
    # predicted occupancy
    y[o,] <- rbinom(species,1,zrep[,site[o],closure[o]]*p[o,])
    
  }
  colnames(y) <- colnames(occ) 

  vis_rep[[i]] <-data.frame(cbind(obs_count= rowSums(y),vis[c("site_5km","TP")],rep=i))
  sum_rep[[i]] <-data.frame(cbind(y_sum=sum(y),rep=i))
  spe_rep[[i]] <-data.frame(cbind(obs_count= colSums(y),rep=i,snames=names(colSums(y))))
  
  saveRDS(y,file=paste0("predicted_observed_occupancy/",data,"y_rep_",i,"_",mod,".rds"))
  

  # could look at transitions states (AHM 2 p243)?
  s_x_e[i] <- sum((colSums(y) - colSums(p))^2/(colSums(p)+e)) # chi square discrepancy over species
  s_x_o[i] <- sum((colSums(occ)  - colSums(p))^2/ (colSums(p)+e))
  v_x_e[i] <- sum((rowSums(y) - rowSums(p))^2/(rowSums(p)+e)) # chi square discrepancy over visits
  v_x_o[i] <- sum((rowSums(occ)  - rowSums(p))^2/ (rowSums(p)+e))
  
  
}

##############################
# look at predictions in space and time
# get site id
site_id <- read.csv("site_id.csv")[2]
colnames(site_id) <- "site_5km"


# get visit info
visit <- jas_data[[2]]
visit$Year <- visit$Year+1993
closure_period <- read.csv("Occ_workflow_V2/date_preparation_1/clsr_per.csv",header=T)[-1]; colnames(closure_period)[2] <- "Year"
closure_period <- cbind(closure_period,Year.x = rep(seq(1994,2019,2),each=2)) 

# set up output data frame
y <- readRDS(paste0("predicted_observed_occupancy/",data,"y_rep_",1,"_",mod,".rds"))

# pivot into longer data set
visit <- left_join(x=visit,y=closure_period,by="Year")
y_long_rep <- cbind(y,visit[c("site_5km","TP","Year.x")]) %>% 
              pivot_longer(colnames(y),names_to = "species") 



# summarize by site and species
# year
y_year_total <- y_long_rep %>% group_by(site_5km,Year.x) %>% summarise(total_y = sum(value))
nr <- nrow(y_year_total)
y_year_out <- data.frame(matrix(nrow=nr,ncol=nrep+2))
y_year_out[,1:3] <- y_year_total[,1:3]

# site
y_site_total <- y_long_rep %>% group_by(site_5km) %>% summarise(total_y = sum(value))
nr <- nrow(y_site_total)
y_site_out <- data.frame(matrix(nrow=nr,ncol=nrep+1))
y_site_out[,1:2] <- y_site_total[,1:2]

# repeat for all repetitions
for(i in 2:nrep){
  # set up output data frame
  y <- readRDS(paste0("predicted_observed_occupancy/",data,"y_rep_",i,"_",mod,".rds"))
  
  # pivot into longer data set
  y_long_rep <- cbind(y,visit[c("site_5km","TP","Year.x")]) %>% pivot_longer(colnames(y),names_to = "species") 
  
  # summarize by site and species
  y_year_out[i+2] <- (y_long_rep %>% group_by(site_5km,Year.x) %>% summarise(total_y = sum(value)))["total_y"]
  y_site_out[i+1] <- (y_long_rep %>% group_by(site_5km) %>% summarise(total_y = sum(value)))["total_y"]
}

# summarize the observed values
obs_y_long <- cbind(occ,visit[c("site_5km","TP","Year.x")]) %>% pivot_longer(colnames(occ),names_to = "species") 
obs_year_total <- obs_y_long %>% group_by(site_5km,Year.x) %>% summarise(total_y = sum(value))
obs_year_total$Year.x <- as.character(obs_year_total$Year.x)
obs_site_total <- obs_y_long %>% group_by(site_5km) %>% summarise(total_y = sum(value)) 
 

# set up with a covariate
# get site id
site_id <- read.csv("site_id.csv")[2]
colnames(site_id) <- "site_5km"

# Check NA values!
# compared predicted with observed distributions
compared_distributions <- function(cov.frame,x.lab,y.lab,bin.size,year){
  
  if(year){
    obs_cov_join <- obs_year_total %>% left_join(x=.,y= cov.frame,by = c("site_5km","Year.x"))
    obs_cov_join <- obs_cov_join[obs_cov_join$Year.x!=2018,]
    obs_cov_tally <- rep(obs_cov_join$value,times=obs_cov_join$total_y)
  } else {
    obs_cov_join <- obs_site_total %>% left_join(x=.,y= cov.frame,by = "site_5km")
    obs_cov_tally <- rep(obs_cov_join$use,obs_cov_join$total_y)
  }

  # create observed histogram along the covariate
  obs_histogram <- ggplot() + geom_histogram(aes(x= obs_cov_tally),color="black",alpha = 0.5,bins = bin.size)
  hist_data <- ggplot_build( obs_histogram)$data[[1]]
  
  # join replicated data with covariate information
  
  if(year){
    colnames(y_year_out)[1:2] <-  c("site_5km","Year.x")
    y_year_out$Year.x <- as.character(y_year_out$Year.x)
    rep_cov_join <- y_year_out %>% left_join(x=.,y=cov.frame,by =  c("site_5km","Year.x"))
  } else {
    colnames(y_site_out)[1] <- "site_5km"
    rep_cov_join <- y_site_out %>% left_join(x=.,y=cov.frame,by =  c("site_5km"))
  }
  rep_count <- data.frame(matrix(ncol = nrep,nrow = bin.size)) 
  
  for(i in 1:200){
    
    # replicate observations by estimated count
    if(year){
      rep_cov_join <-  rep_cov_join[rep_cov_join$Year.x!=2018,] 
      rep_cov <- rep(rep_cov_join$value,rep_cov_join[,paste0("X",i+2)])
    } else {
      rep_cov <- rep(rep_cov_join$use,rep_cov_join[,paste0("X",i+1)])
    } 
    
    # use histogram
    rep_histogram <- ggplot() + geom_histogram(aes(x=rep_cov),bins=bin.size)
    hist_data_rep <- ggplot_build(rep_histogram)$data[[1]]
    
    # check the bins match from the observed and replicated data
    if(all(hist_data$x!=hist_data_rep$x)){stop()}
    rep_count[,i] <- hist_data_rep$y
  }  
  
  # bind with x values
  # comparative histogram
  hist_summary_data <- data.frame(summary_f(rep_count),xval=hist_data$x)
  out <- obs_histogram + geom_point(aes(y=hist_summary_data$mean.x,x=hist_summary_data$xval),color="blue") +
                  geom_line(aes(y=hist_summary_data$mean.x,x=hist_summary_data$xval),color="blue") + 
                  geom_errorbar(aes(ymin=hist_summary_data$ql,ymax=hist_summary_data$qu,x=hist_summary_data$xval),color="blue") +
                  xlab(x.lab) + ylab(y.lab) +
                  theme( panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         panel.background = element_blank(), 
                         axis.text.y = element_text(size=10), 
                         axis.text.x= element_text(size=10), 
                         axis.title.x= element_text(size=11), 
                         axis.title.y =element_text(size=11), 
                         axis.line = element_line(colour = "black")) 
  print(out)
  return(list(plot=out,observed_data=hist_data,predicted_data=rep_count))
}

# aggregated observed values vs predicted values
covs_observed <- readRDS("covars_wide_bees.rds")

# mean temperature
cov.frame  <- covs_observed$mean_temp 
colnames(cov.frame)[2] <- "use"
colnames(cov.frame)[1] <- "site_5km"  
mean_temp_pp <- compared_distributions(cov.frame,x.lab = "Spatial Temperature",y.lab = "Observation Count",bin.size = 40,year=F)

# mean risk quotient
cov.frame = covs_observed$RQsum_M
colnames(cov.frame)[2] <- "use"
colnames(cov.frame)[1] <- "site_5km"   
mean_rqm_pp <- compared_distributions(cov.frame,x.lab = "Spatial Risk Quotient",y.lab = "Observation Count",bin.size = 40,year=F)

# mean agriculture
cov.frame =  covs_observed$agri
colnames(cov.frame)[2] <- "use"
colnames(cov.frame)[1] <- "site_5km"   
mean_agri_pp <- compared_distributions(cov.frame,x.lab = "Agriculture landcover",y.lab = "Observation Count",bin.size = 40,year=F)

# mean agriculture
cov.frame =  covs_observed$semi
colnames(cov.frame)[2] <- "use"
colnames(cov.frame)[1] <- "site_5km"   
mean_semi_pp <- compared_distributions(cov.frame,x.lab = "Semi-natural landcover",y.lab = "Observation Count",bin.size = 40,year=F)

# temporal temperature
cov.frame =  pivot_longer(covs_observed$temp_anom,colnames(covs_observed$temp_anom[-1]),names_to = "Year.x")
colnames(cov.frame)[1] <- "site_5km"   
mean_tmpa_pp <- compared_distributions(cov.frame,x.lab = "Temporal Temperature",y.lab = "Observation Count",bin.size = 40,year=T)

# temporal risk quotient
cov.frame =  pivot_longer(covs_observed$temp_anom,colnames(covs_observed$RQsum_A[-1]),names_to = "Year.x")
colnames(cov.frame)[1] <- "site_5km"   
mean_rqa_pp <- compared_distributions(cov.frame,x.lab = "Temporal Risk Quotient",y.lab = "Observation Count",bin.size = 40,year=T)

cov_pp_plot <- ggarrange(mean_temp_pp$plot,
          mean_tmpa_pp$plot,
          mean_agri_pp$plot,
          mean_semi_pp$plot,
          mean_rqm_pp$plot,
          mean_rqa_pp$plot,ncol=2,nrow=3)

ggsave(cov_pp_plot,"cov_pp_plot.png")


# use a histogram
##############################
# look at predicted values vs. real values
species_p <- mean(s_x_e>s_x_o)
visit_p <- mean(v_x_e>v_x_o)

xys  <- seq(min(c(s_x_e,s_x_o))-5, max(c(s_x_e,s_x_o))+5 )
xyv  <- seq(min(c(v_x_e,v_x_o))-5, max(c(v_x_e,v_x_o))+5 )

mins <- min(c(s_x_e,s_x_o))

visit_pp_plot <- ggplot() + geom_point(aes(v_x_o,v_x_e )) + geom_line(aes(xyv,xyv)) + ggtitle(paste0(mod," ","Visit B p-value = ",visit_p ))
species_pp_plot <- ggplot() + geom_point(aes(s_x_o,s_x_e )) + geom_line(aes(xys,xys)) + ggtitle(paste0("Species B p-value = ",species_p))

ggsave(paste0("model_fit/", data,"_",mod,".png") ,ggarrange(visit_pp_plot,species_pp_plot),width = 7,height = 5)

# put lists into date frames
vis_rep. <- do.call(rbind, vis_rep)
sum_rep. <- do.call(rbind, sum_rep)
spe_rep. <- do.call(rbind, spe_rep)

spe_rep.$obs_count <- as.numeric(spe_rep.$obs_count )


act <-colSums( occ)


# plot all species fit
nm <- unique(spe_rep.$snames)
names(nm) <- unique(spe_rep.$snames)
splots <- lapply(nm,function(x){
  
  ggplot() + geom_histogram(data=spe_rep.[spe_rep.$snames==x,] ,aes(x=obs_count))+
    geom_vline(aes(xintercept=act[names(act)==x]),color="blue")+ggtitle(paste(x))
})

indplots <-ggpubr::ggarrange(plotlist=splots,ncol=4,nrow=4)

for(i in 1:length(indplots)){
  
ggsave(paste0("model_fit/",mod,"_",data,"_species_fit_",i,".png"), indplots[[i]],width=10,height=8)
  
}
   
# overall model fit
sum_p <- ggplot() + geom_histogram(data= sum_rep.
                                   ,aes(x=y_sum),bins =40)+
  geom_vline(aes(xintercept=sum(rowSums(occ))),color="blue")+ggtitle(paste0('Total Count Estimates \nModel = ',data," ",mod))

ggsave(paste0("model_fit/",mod,"_",data,"_all_fit.png"), sum_p,width=4,height=4 )

###############################################################
##################### Run simulations #########################
###############################################################
# ALL MEAN EFFECT
out <- occ_output[[2]]
# iterations used
iterations. <- 2000

# beta coefficients
beta <- out$sims.list$mu.beta[1:iterations.,]



# gamma
gamma <- inv_logit_scaled( out$sims.list$mu.gamma[1:iterations.] )


# initial occupancy
init <- out$sims.list$init.occ

# alpha.phi - intercept
alpha.phi <- out$sims.list$mu.alpha.phi[1:iterations.]
 

# RISK QUOTIENT SPATIAL EFFECT
# covariates
covs <- c(base) 
covs1 <- c(base)

# cov 1
# min vlue
covs$RQA[,1:12] <- 0 # set temporal effect to zero

# cov 2
# min vlue
min(covs$RQA)
covs1$RQA[,1:12] <- 0 # set temporal effect to zero

min(covs$RQM)
covs1$RQM[,1:12] <-  -5.855139 # centered value for zero

# simulation
psi1 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                beta = beta,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs),
                species=F)

psi2 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                beta = beta,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs1),
                species=F
                )

# summarize all sites by year
m1 <- summ(psi1,year=seq(1994,2018,2))
m1$name <- "With pesticide"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Without pesticide"

m <- rbind(m1,m2)

# occupancy by year
ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name),size=1)+
         geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")

# difference risk quotient scenarios
diffp <-  psi1-psi2

 
diffsum <- summ(diffp,year = seq(1994,2018,2)) 

# difference plot
ggplot(data=diffsum)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")


# site level changes between the first and the last year
sitep <- summs(diffp,year = 13)

# read in the sight id grid cells
site_id <- read.csv("site_id.csv")[2]
site_changes <- cbind(sitep,site_id)
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
sitech <- left_join(site_changes,coord )

# rq effect map
plots_map <- list(
  ggplot() +
    geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(90000, 660000) +ylim(0,660000)+ 
    geom_tile(data = sitech, 
              aes(x = E, y = N, fill =mean))+
    scale_fill_continuous(type = "viridis", name = "",direction=-1)+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                           axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                           axis.title.x=element_blank(),
                                                                           legend.title=element_blank(),
                                                                           axis.title.y=element_blank(),
                                                                           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                           panel.grid.minor=element_blank(),plot.background=element_blank())+
    ggtitle("Mean change in occupancy"),
  
  ggplot() +
    geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
    ylim(0, 700000)  + geom_tile(data = sitech, 
                                 aes(x = E, y = N, fill =lb90))+
    scale_fill_continuous(type = "viridis", name = "")+ theme( axis.text.y = element_text(size=10),
                                                               axis.text.x=  element_text(size=10),
                                                               axis.title.x= element_text(size=11), 
                                                               axis.title.y =element_text(size=11),panel.grid.major = element_blank(), 
                                                               panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(),
    )+
    ggtitle("Lower 95%CI change in occupancy"),
  
  ggplot() +
    geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
    ylim(0, 700000)  + geom_tile(data = sitech, 
                                 aes(x = E, y = N, fill =ub90))+
    scale_fill_continuous(type = "viridis", name = "")+ theme( axis.text.y = element_text(size=10),
                                                               axis.text.x=  element_text(size=10),
                                                               axis.title.x= element_text(size=11), 
                                                               axis.title.y =element_text(size=11),panel.grid.major = element_blank(), 
                                                               panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(),
    )+
    ggtitle("Upper 95%CI change in occupancy")
)
ggpubr::ggarrange(plotlist = plots_map,nrow=1)
ggsave("spatial_plot_mean.png",dpi=2000,plot= plots_map[[1]],width = 4.5,height=4.5)
ggsave("spatial_plot.png",dpi=2000,plot= plots_map[[2]],width = 4.5,height=4.5)




################################################################################
# TEMPORAL RISK QUOTIENT EFFECT
# covariates
covs <- c(base)  # original covariates
covs1 <- c(base) # covariates to be changed


# covariates
rqtemp <-  readRDS("zero_app.rds") # standardized value for zero for each covariate

covs$RQA[,1:12] <- rqtemp$actual[-1] 
covs1$RQA[,1:12] <-  rqtemp$zero_application[-1]


# Pesticide Simulation 
# parameters for simulation
psi1 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                beta =  beta,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs),
                species=F)

sim.pop(psi.start = init,
        gamma = gamma,
        beta = beta,
        alpha = alpha,
        psi =  array(dim=c(nrep,nspecies,sites,time)),
        cov.array = f_array(covs),
        species=T,
        region = region.psi,
        region.cov = jas_data$region
)
psi2 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                beta =  beta,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs1),
                species=F)


# different between simulations
diffp<-  psi1-psi2


# summarize occupancies
m1 <- summ(psi1,year=seq(1994,2018,2))
m1$name <- "With pesticide"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Without pesticide"


m3 <- rbind(m1,m2)
rqt_summ <- summ(diffp,year = seq(1994,2018,2)) 

# rq temporal effect through time
ggplot(data=rqt_summ)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")



###############################################################
# site level changes between for all years
plots <- list()
years <- seq(1996,2018,2)
for(i in c(2,12,13)){
  sitep <- summs(diffp,year = i)
  
  site_id <- read.csv("site_id.csv")[2]
  sitech <- cbind(sitep,site_id)
  coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
  sitech <- left_join(sitech,coord )
  sitech$mean[sitech$ub95>=0] <- NA
  plots[[i-1]] <- ggplot() +
    geom_path(data = UK$britain, aes(x = long, y = lat, group = group))+ xlim(90000, 660000) +ylim(0,660000)+
    geom_tile(data = sitech, 
              aes(x = E, y = N, fill =mean))+
    scale_fill_continuous(type = "viridis", name = paste(years[i-1]),na.value="grey",limits=c(-0.095,0))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          legend.text = element_text(size=12),
          legend.title= element_text(size=12),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
}

# temporal maps
tmaps<- ggarrange(plotlist=plots,nrow=2,ncol=3) # map of all temporal changes
map<- ggarrange(plotlist=list(plots[[1]],plots[[11]],plots[[12]]) ,nrow=1,ncol=3)
ggsave("temp_plot1.png",plot=  tmaps$`1`,width = 10,height=7)
ggsave("temp_plot2.png",plot= tmaps$`2`,width = 10,height=7)
ggsave("temporal_rq.png",plot= map,width = 16,height=7,dpi = 2000)

  
############################################################################
##################### Persistence (%) Marginal Effects #####################
############################################################################
covsr <- readRDS("covars_wide.rds")

# temperature spatial marginal effect ######################################
temps <- covsr$temp_anom[-1]+covsr$mean_temp[-1] # non standardized temperatures

me <-mean(as.matrix( covsr$mean_temp[-1])) # mean
mn<- min(covsr$mean_temp[-1])-me # centered minimum temperature
mx<- max(covsr$mean_temp[-1])-me # centered minimum temperature

# intercept
alpha = alpha.phi

# values used to predict
pout <- seq(mn,mx, length.out =100)

# predict persistence
phi <-  apply( t(matrix(ncol=iterations.,nrow=100,(pout))) *beta[,1] ,2,function(x){inv_logit_scaled( x+alpha)})

sphi <- summp(phi)

# plot persistence
tmp <- ggplot(data=sphi)+geom_line(aes(x=pout+me,y=mean*100),size=1)+
  geom_ribbon(aes(x=pout+me,ymin=lb*100,ymax=ub*100),alpha=0.2)+
  geom_vline(xintercept =me,color="black",linetype="dashed") +
  
  theme(        axis.text.y = element_text(size=10),
                axis.text.x=  element_text(size=10),
                axis.title.x= element_text(size=11), 
                axis.title.y =element_text(size=11),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Mean Temperature (C)")

ggsave("tmpr_eff.png",plot= tmp,width = 6,height=3)

# temperature temporal marginal effect #########################################
temps <- covsr$temp_anom[-1] # non standardized temperatures

me <- mean(as.matrix(temps)) # mean
min <- min(temps)
max <- max(temps)

mn <- min # centered minimum temperature
mx <- max 

# intercept
alpha = alpha.phi

# values used to predict
pout <- seq(mn,mx, length.out =100)

# predict persistence
phi <-  apply( t(matrix(ncol=iterations.,nrow=100,(pout))) *beta[,2] ,2,function(x){inv_logit_scaled( x+alpha)})

sphi <- summp(phi)

# plot persistence
tmp <- ggplot(data=sphi)+geom_line(aes(x=pout+me,y=mean*100),size=1)+
  geom_ribbon(aes(x=pout+me,ymin=lb*100,ymax=ub*100),alpha=0.2)+
  geom_vline(xintercept =me,color="black",linetype="dashed") +
  
  theme(        axis.text.y = element_text(size=10),
                axis.text.x=  element_text(size=10),
                axis.title.x= element_text(size=11), 
                axis.title.y =element_text(size=11),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Mean Temperature (C)")

ggsave("tmpr_eff.png",plot= tmp,width = 6,height=3)

# risk temporal marginal effect ##############################################
rq <- covsr$RQsum_A[-1]+covsr$RQsum_M[-1] # non standardized risk quotient

me=mean(as.matrix(rq)) # mean
mx=max(rq)-mean(as.matrix(rq))  # centred maximum
mn=min(rq)-mean(as.matrix(rq)) # centred minimum


per<-(((seq(mn,mx, length.out =200)+me)-me)/me)*100 # anomalies expressed as a percentage

# calculate persistence
phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,5] ,2,function(x){inv_logit_scaled( x+alpha)})

# summarize persistence
spout <- summp(phi)
nf<- quantile(((rq/covsr$RQsum_M[-1])/covsr$RQsum_M[-1])*100,0.95,na.rm=T) # upper 95CI

me=mean(as.matrix(rq))

rqap<-ggplot(data=spout)+geom_line(aes(x=per,y=mean*100),size=1)+
  geom_ribbon(aes(x=per,ymin=lb*100,ymax=ub*100),alpha=0.2)+
  geom_vline(xintercept =nf,color="red",linetype="dashed") +
  geom_vline(xintercept =-100,color="red",linetype="dashed") +
  geom_vline(xintercept =0,color="black",linetype="dashed") +
  theme(    axis.text.y = element_text(size=10),
            axis.text.x=  element_text(size=10),
            axis.title.x= element_text(size=11), 
            axis.title.y =element_text(size=11),
            
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Percentage change in risk quotient at a site (%)")

ggsave("rqap.png",plot= rqap,width = 6,height=5)

# semi natural marginal effect ##############################################
semi <- covsr$semi[-1]
mx=max(semi)-mean(as.matrix(semi))
mn=min(semi)-mean(as.matrix(semi))
max <-max(semi)
min <-min(semi)
phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,3] ,2,function(x){inv_logit_scaled( x+alpha)})

sephi <- summp(phi)


me=mean(as.matrix(semi))

semip <- ggplot(data=sephi)+geom_line(aes(x=seq(min,max, length.out =200),y=mean*100),size=1)+
  geom_ribbon(aes(x=seq(min,max, length.out =200),ymin=lb*100,ymax=ub*100),alpha=0.2)+
  geom_vline(xintercept =me,color="black",linetype="dashed") +
  
  theme(        axis.text.y = element_text(size=10),
                axis.text.x=  element_text(size=10),
                axis.title.x= element_text(size=11), 
                axis.title.y =element_text(size=11),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Percentage semi-natural land cover")
semip
ggsave("semi_eff.png",plot= semip,width = 6,height=3)


# agriculture marginal effect #######################################

agri <- covsr$agri[-1]
mx=max(agri)-mean(as.matrix(agri))
mn=min(agri)-mean(as.matrix(agri))
max=max(agri)
min=min(agri)
me=mean(as.matrix(agri))
agphi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,3] ,2,function(x){inv_logit_scaled( x+alpha)})
agphi <- summp(agphi)

agrip <- ggplot(data=agphi)+geom_line(aes(x=seq(min,max, length.out =200),y=mean*100),size=1)+
  geom_ribbon(aes(x=seq(min,max, length.out =200),ymin=lb*100,ymax=ub*100),alpha=0.2)+
  geom_vline(xintercept =me,color="black",linetype="dashed") +
  
  theme(        axis.text.y = element_text(size=10),
                axis.text.x=  element_text(size=10),
                axis.title.x= element_text(size=11), 
                axis.title.y =element_text(size=11),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Agricultural land cover")
agrip

######################
rq <- covsr$agri[-1]
mx=max(rq)-mean(as.matrix(rq))
mn=min(rq)-mean(as.matrix(rq))
max=max(rq)
min=min(rq)
phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,4] ,2,function(x){inv_logit_scaled( x+alpha)})

b1.1 <- summp(phi)
b2 <- summp(data.frame( inv_logit_scaled(alpha)))

b1$Landcover <- "Semi-natural"
b1$x <- seq(min,max, length.out =200)
b1.1$Landcover <- "Agriculture"
b1.1$x <- seq(min,max, length.out =200)


me=mean(as.matrix(rq))


b1 <- rbind(b1,b1.1)
b1$Landcover <-as.factor(b1$Landcover)

ggplot(data=b1)+geom_line(aes(x=x,y=mean , color=Landcover )  ,size=1)+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,fill=Landcover ),alpha=0.2)+
  theme(        axis.text.y = element_text(size=10),
                axis.text.x=  element_text(size=10),
                axis.title.x= element_text(size=11), 
                axis.title.y =element_text(size=11),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
  )+ylab("Persistence (%)")+xlab("Percentage land cover")
ggsave("semi_eff.png",plot= semi,width = 6,height=5)


