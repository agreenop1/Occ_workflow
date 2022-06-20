################################################################################
################################################################################
##################### simulate different scenarios #############################
library(brms)
library(BRCmap)
library(ggmcmc)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
source("Occ_workflow_V2/outputs_3/combine_chain_4.2.R")
UK <-  readRDS("UK_map.rds")
# read in parameters
para <- c("mu.beta","gamma","beta","init.occ","mu.alpha.phi","mu.gamma","alpha.phi",
          "dtype1.p","dtype2.p","dtype3.p","alpha.p")

# start of file name
file="3_bee" 

# read in data with cover it information
data="bees"


# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))

# the number of sites and number of time periods
time = 13
sites = 2278
species = ncol(jas_data[[1]][-1])
# number of iterations in each file
its <- c(500)

# file path where the files are located
file.path = paste0("jas_out/",file ,"_C.CID_ID_UID.rdata")

# iterations used
iterations. <- 1800


# covariates
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)


# compile the chains from the model
#out <- comb_daisy(parameters=para,
#                  iter.index=33:33,chain.index=1:3,summary=T,file.path=file.path,by.it=500,
#                  it.used=its,
#                  iterations=1000,verbose=T)
#
occ_output <- readRDS("jas_out/summary_ps.rds")
out <- occ_output[[2]]



############################## Array function ################
# create array for covs
f_array <- function(x){
             row = nrow(x[[1]]); col = ncol(x[[1]]); nmat = length(covs)
             cov.array <- array(dim=c(row,col,nmat))
             
             for(n in 1:nmat){
               cov.array[,,n]  <- as.matrix(x[[n]])}
cov.array             
}
##############################################################
################# Simulation function ########################
##############################################################
sim.pop <-  function(psi,psi.start,alpha,gamma,beta,cov.array,species){

if(species==T){

nspecies <- ncol(gamma)

for(s in 1:nspecies){  
  
  # starting values
  psi[,s,1:sites,1] <- psi.start[,s]
  
  # calculate phi
  for(i in 1:sites){
    for(t in 2:time){
    
      
      phi <- alpha[,s] + cov.array[i,t-1,1]*beta[,1,s] +
                     cov.array[i,t-1,2]*beta[,2,s] +
                     cov.array[i,t-1,3]*beta[,3,s] +
                     cov.array[i,t-1,4]*beta[,4,s] +
                     cov.array[i,t-1,5]*beta[,5,s] +
                     cov.array[i,t-1,6]*beta[,6,s] 
      
      psi[,s,i,t] <-  psi[,s,i,t-1]* inv_logit_scaled(phi) + (1- psi[,s,i,t-1])*gamma[,s]
      
    }
  }
}  
psi

}else{
  # starting values
  psi[,1:sites,1] <- psi.start
  
  # calculate phi
  for(i in 1:sites){
    for(t in 2:time){
      
      
      phi <- alpha + cov.array[i,t-1,1]*beta[,1] +
        cov.array[i,t-1,2]*beta[,2] +
        cov.array[i,t-1,3]*beta[,3] +
        cov.array[i,t-1,4]*beta[,4] +
        cov.array[i,t-1,5]*beta[,5] +
        cov.array[i,t-1,6]*beta[,6] 
      
      psi[,i,t] <-  psi[,i,t-1]* inv_logit_scaled(phi) + (1- psi[,i,t-1])*gamma
    }
  }
  psi 
  
}
}

# summarize output
summ <- function(x,year){
  m1 <- apply(x,c(3),mean)
  m2 <- apply(x,c(3),quantile,probs = 0.025)
  m3 <- apply(x,c(3),quantile,probs = 0.975)
  
  trend <- data.frame(mean=m1,lb=m2,ub=m3)
  trend$year <- year 
  trend}



#################################################################
############# PP Checks #########################################
#################################################################
# State Model
# aggregate observations where a species was observed at a site

occ <- jas_data[[1]][-1]
vis <- jas_data[[2]]
occdat <- cbind(occ,vis) 



# simulation
covs <- c(base) # covariate list


#  multiple species check
nrep <- 50 #dim(out$sims.list$beta)[1]
popocc <- array(dim=c(nrep,species,sites,time)) # occupancy estimates
beta <- out$sims.list$beta[1:nrep,,] # beta coefficient
gamma <-  out$sims.list$gamma[1:nrep,] # colonization
init <-  out$sims.list$init.occ[1:nrep,] # initial occupancy
alpha <- out$sims.list$alpha.phi[1:nrep,] # intercepts


# estimates of population occupancy for species
psi <- sim.pop(psi.start = init,
                gamma = gamma,
                beta = beta,
                alpha = alpha,
                psi =  array(dim=c(nrep,species,sites,time)),
                cov.array = f_array(covs),
                species=T)

# predict occupancy state
popbin <- apply(psi,c(2,3,4),function(x){rbinom(length(x),1,x)})

# observed occupancy
z <- jas_data[[3]]

# observation model
out <- occ_output[[3]]
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
y <- array(dim=c(nrep,nobs,species)) # predicted observations
p <- array(dim=c(nrep,nobs,species)) # probability occupancy


# observation model
for(i in 1:nrep){
  for(o in 1:nobs){
 
    # observation probability  
    p[i,o,] <- inv_logit_scaled( alpha.p[i,closure[o]] + d1[i,] + d2[i,]*SHORT[o] +
      d3[i,]*LONG[o])
    
    # predicted occupancy
    y[i,o,] <- rbinom(94,1,popbin[i,,site[o],closure[o]]*p[i,o,])
    
  }
  colnames(y[i,,]) <- colnames(occ)
}

# aggregate all observations seen at a site
yrep <- cbind(obs_count= rowSums(y[1,,]),vis[c("site_5km","TP")],rep=1)

for(i in 1:nrep){
        
      yreps <- cbind(obs_count= rowSums(y[i,,]),vis[c("site_5km","TP")],rep=i)
      
      yrep <- rbind(yrep,yreps)
                       
}

# look at predicted values vs. real values
ggplot() + geom_density(data=yrep,aes(x=obs_count,group=as.factor(rep)),adjust=2.5)+
  geom_density(aes(x=rowSums(occ)),color="blue", adjust=2.5)


# state model only
# get where species were observed at a site
check <- cbind(occ,vis[c("site_5km","TP")]) %>% group_by(site_5km,TP) %>% summarise_all(sum) # observed occupancy
check[3:96][check[3:96]>0] <- 1
sr_site <- cbind(total=rowSums( check[3:96]),check[c("site_5km","TP")]) # observed occupancy sr

# summarise predicted occupancies status
pop1 <- popbin[1,,,1] # get occupancy status
zi <- z[,,1] # observed occupancies
zi[is.na(zi)] <- 0 # set any na to zero
pop1[zi==0] <- 0 # make sure we only include observations where the species has been observed

occres1 <- data.frame(sr=colSums(pop1),sr_diff=colSums(pop1)-colSums(zi),rep=1,tp=1)

# repeat for reps 
for(i in 1:nrep){
  for(t in 1:time){
    
   pop1 <- popbin[i,,,t]
   zi <- z[,,t]
   zi[is.na(zi)] <- 0
   pop1[zi==0] <- 0
   occres <- data.frame(sr=colSums(pop1),sr_diff=colSums(pop1)-colSums(zi),rep=i,tp=t)
   occres1 <- rbind(occres1,occres)
   
  }
}

# plots of model predictive ability
ggplot() + geom_density(data=occres1,aes(x=log(sr),group=as.factor(rep)))+geom_density(aes(x=log(sr_site$total)),color="blue")
ggplot() + geom_density(data=occres1,aes(x=sr_diff,group=as.factor(rep)))


###############################################################
##################### Diagnostic Plots ########################
###############################################################
source("output_functions.R")
species <- colnames(jas_data[[1]][-1])
nspecies <- length(species)
out <- occ_output[[2]]

# names of covariates
covs <- c("Temperature spatial","Temperature temporal",
          "Semi Natural Landcover","Agricultural Landcover"
          ,"Risk Quotient temporal","Risk Quotient spatial"
)

n.cov <- length(covs)

# parameters to check
parameters=c("mu.beta")

# output <- summary_chains(out,comb.chain = T,keep.samples = F)
options(max.print=10000);print(out)

# model 
samples = out$samples
occ.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,] # check parameter convergence
rhat. <- out$Rhat # all rhat

samp.df <- ggmcmc::ggs(samples) # all samples

# par.p plots all parameters above - mainly diagnostics -  should detect each type of parameter okay
# parameters need to be a vector of names
par.p <-sapply(parameters,output_plots,USE.NAMES=T)
par.p


# save plots
wid = 6
hei = 3

par.p$mu.beta

# plots of chains, density and parameter estimate 

m.beta <- ggpubr::ggarrange(plotlist = par.p$mu.beta.mu.beta,nrow=1,ncol=2)
ggsave(paste0("effect_plot/Main Effects.1",'_',file,".png"),m.beta[[1]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.2",'_',file,".png"),m.beta[[2]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.3",'_',file,".png"),m.beta[[3]],width = wid,height = hei)
m.beta

# individual species plots 
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



###############################################################
##################### Run simulations #########################
###############################################################
#   
# out latent occupancy
beta <- out$sims.list$mu.beta
ecolMeans(beta)
ncol(s)
# gamma
gamma <- inv_logit_scaled( out$sims.list$mu.gamma )
summary(mean(gamma) )

# initial occupancy
init <- out$sims.list$init.occ
inv_logit_scaled (mean(colMeans(logit_scaled( init))))

# alpha.phi - intercept
alpha.phi <- out$sims.list$mu.alpha.phi
summary(mean(alpha.phi)) 

# SPATIAL
# Pesticide Simulation 
# parameters for simulation
covs <- c(base) 
covs1 <- c(base)

# cov 1
# min vlue
min(covs$RQA)
covs$RQA[,1:12] <- 0

min(covs$RQM)


# cov 2
# min vlue
min(covs$RQA)
covs1$RQA[,1:12] <- 0

min(covs$RQM)
covs1$RQM[,1:12] <-  -5.855139

# simulation
psi1 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs))

psi2 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs1))


m1 <- summ(psi1,year=seq(1994,2018,2))
m1$name <- "With pesticide"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Without pesticide"

m <- rbind(m1,m2)

library(ggplot2)
ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name),size=1)+
         geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")

diffp <-  psi1-psi2

 
 m3 <- summ(diffp,year = seq(1994,2018,2)) 


ggplot(data=m3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")



###############################################################
# Change Maps 
# summarize output
summs <- function(x,year){
  m1 <-round( apply(x,c(2,3),mean),3)
  m2 <-round( apply(x,c(2,3),quantile,probs = 0.025),3)
  m3 <-round( apply(x,c(2,3),quantile,probs = 0.975),3)
  m4 <-round( apply(x,c(2,3),quantile,probs = 0.05),3)
  m5 <-round( apply(x,c(2,3),quantile,probs = 0.95),3)
  m6 <-round( apply(x,c(2,3),sd),3)
if(is.null(year)){
list(mean = m1,lb = m2, ub = m3)}else{
  data.frame(mean = m1[,year],lb95 = m2[,year], ub95 = m3[,year],
             lb90 = m4[,year], ub90 = m5[,year],sd=m6[,year])
}
}



# site level changes between the first and the last year
sitep <- summs(diffp,year = 13)

site_id <- read.csv("site_id.csv")[2]
sitech <- cbind(sitep,site_id)
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
sitech <- left_join(sitech,coord )

plots_map <- list(
ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = sitech, 
                               aes(x = E, y = N, fill =mean))+
  scale_fill_continuous(type = "viridis", name = "")+ theme(panel.grid.major = element_blank(), 
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.background = element_blank())+
                                                                       ggtitle("Mean change in occupancy"),

ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = sitech, 
                               aes(x = E, y = N, fill =sd))+
  scale_fill_continuous(type = "viridis", name = "")+ theme(panel.grid.major = element_blank(), 
                                                                                    panel.grid.minor = element_blank(),
                                                                                    panel.background = element_blank())+
  ggtitle("SD change in occupancy")
)

ggsave("spatial_plot.png",plot= ggpubr::ggarrange(plotlist = plots_map),width = 10,height=7)



###############################################################
##################### Run simulations #########################
###############################################################
# TEMPORAL
# covariates
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)
covs <- c(base) 
covs1 <- c(base)

# cov 1
# min vlue
all(covs$RQA>=covs1$RQA)

# cov 2
# min vlue
min(covs$RQA)
min(covs$RQM)
z0 <-  readRDS("zero_app.rds")
z0$z0>z0$z1
covs$RQA[,1:12] <- z0$z0[-1]
covs1$RQA[,1:12] <-  z0$z1[-1]


# Pesticide Simulation 
# parameters for simulation
psi1 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs))

psi2 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs1))

diffp<-  psi1-psi2
all0 <- distinct(data.frame( apply(diffp,c(1,2),function(x) all(x==0))))
l <- !unlist( all0[1,])
nr <-  nrow(diffp[1,l,])
difference <- array(dim=c(1800,nr,13))
                    
for (i in 1:1800){
difference[i,,]  <- diffp[i,nr,] 
}


m1 <- summ(psi1,year=seq(1994,2018,2))
m1$name <- "With pesticide"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Without pesticide"

m <- rbind(m1,m2)

library(ggplot2)
ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name),size=1)+
  geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")

slope <- function(x){
  l <- dim(x)[1]
  time <- dim(x)[3]
  
  out <- array(dim=c(l,dim(x)[2],time-1))
  
  for(i in 1:l){
    
    out[i,,]  <- (x[i,,-1]-x[i,,1:12])/2
    
  }
  out
}

s1 <- slope(psi1)
s2 <- slope(psi2)

diff_s <- s1-s2
diffp<-  psi1-psi2
m2 <- summ(diff_s,year = seq(1996,2018,2)) 
m3 <- summ(diffp,year = seq(1994,2018,2)) 

ggplot(data=m2)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in slopes")+xlab("Year")

ggplot(data=m3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")


###############################################################
# site level changes between the first and the last year
plots <- list()

for(i in 1:12){
  sitep <- summs(diff_s,year = i)
  
  site_id <- read.csv("site_id.csv")[2]
  sitech <- cbind(sitep,site_id)
  coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
  sitech <- left_join(sitech,coord )
  sitech. <- sitech[sitech$ub95<0,]
  plots[[i]] <- ggplot() +
    geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
    ylim(0, 700000)  + geom_tile(data = sitech., 
                                 aes(x = E, y = N, fill =mean))+
    scale_fill_continuous(type = "viridis", name = paste("Change in occupancy",i))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

plots


###############################################################
# Climate Simulation 
covs2 <- c(base)
covs2$temp.a <- data.frame(covs2$temp.a)
l1 <- covs2$temp.a>0
covs2$temp.a[l1] <- covs2$temp.a[l1]*0.1
 
# parameters for simulation
psi1 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs))

psi2 <- sim.pop(psi.start = 0.50,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs2))


c1 <- summ(psi1,year=seq(1994,2018,2))
c1$name <- "With pesticide"
c2 <- summ(psi2,year = seq(1994,2018,2)) 
c2$name <- "Without pesticide"

c <- rbind(c1,c2)

library(ggplot2)
ggplot(data=c)+geom_line(aes(x=year,y=mean,color=name),size=1)+
  geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")

slope <- function(x){
  l <- dim(x)[1]
  time <- dim(x)[3]
  
  out <- array(dim=c(l,dim(x)[2],time-1))
  
  for(i in 1:l){
    
    out[i,,]  <- (x[i,,-1]-x[i,,1:12])/2
    
  }
  out
}

s1 <- slope(psi1)
s2 <- slope(psi2)

diff_s <- s1-s2
diffc<-  psi1-psi2
c2 <- summ(diff_s,year = seq(1996,2018,2)) 
c3 <- summ(diffc,year = seq(1994,2018,2)) 

ggplot(data=c2)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in slopes")+xlab("Year")

ggplot(data=c3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")
############################################################################

# site level changes between the first and the last year
sitec <- summs(diffc,year = 13)

site_id <- read.csv("site_id.csv")[2]
sitecc <- cbind(sitec,site_id)
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
sitecc <- left_join(sitecc,coord )

ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = sitecc, 
                               aes(x = E, y = N, fill =mean))+
  scale_fill_continuous(type = "viridis", name = "Change in occupancy")+ theme(panel.grid.major = element_blank(), 
                                                                                     panel.grid.minor = element_blank(),
                                                                                     panel.background = element_blank())





