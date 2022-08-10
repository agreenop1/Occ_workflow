################################################################################
################################################################################
##################### simulate different scenarios #############################
library(brms)
library(ggmcmc)
library(ggplot2)
library(stringr)
library(pROC)
library(dplyr)
library(forcats)
library(caret)
source("output_functions.R")
source("Occ_workflow_V2/outputs_3/combine_chain_4.2.R")

# uk map
UK <-  readRDS("UK_map.rds")
# sre = without year random effect; rre = sum to zero random effect;jd = new year effect not sum zero; new year effect sum zero
# read in data with cover it information
data="hoverflies"




tribe= "crop"
rq = "ld_jdrre"
group = paste0('hoverflies_',rq)
mod =paste0( "all_",tribe)
cat(mod,group,"\n")


hover <-   readRDS("hover_cropSpecies.rds")

#table(hover$tribe)
#genus <- distinct(data.frame(genus=str_split(hover$species," ",simplify = T)[,1],tribe=hover$tribe))
#write.csv(genus,"genus.csv")
unique(hover$tribe)
keep_names <-  readRDS("hover_cropSpecies.rds")

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016","pol",".rds"))


occ <- jas_data[[1]][-1]
# remove id column
occ <- occ[c(keep_names)]
vis <- jas_data[[2]]
occdat <- cbind(occ,vis) 


# the number of sites and number of time periods
time = 13
sites = nrow(jas_data[[4]]$temp_anom)
species = ncol(occ)

# covariates
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)


# output from jags  
occ_output <- readRDS(paste0("jas_out/all_",tribe,"_hoverflies_",rq,'_summary.rds'))
out <- occ_output

cat((nrow(out$sims.list$alpha.phi)*5)/3,"iterations used before thin per chain","\n")
################################################################################
########################### Diagnostic Plots ###################################
################################################################################

# model 
species <- colnames(occ)
nspecies <- length(species)
out <- occ_output

# names of covariates
covs <- c("Temperature spatial","Temperature temporal",
          "Semi Natural Landcover","Agricultural Landcover"
          ,"Risk Quotient temporal","Risk Quotient spatial"
)

n.cov <- length(covs)

# parameters to check
parameters=c("mu.beta","beta","alpha.p")#,
          #   ,
           #  "dtype1.p","dtype2.p","dtype3.p",
           #  "gamma")


# options(max.print=10000);print(out)

samples = out$samples
occ.sum <- data.frame(round(out$summary,3))
con.f <- occ.sum[occ.sum$Rhat>1.05,] # check parameter convergence
rhat <- out$Rhat # all rhat

samp.df <- ggmcmc::ggs(samples) # all samples

# par.p plots all parameters above - mainly diagnostics -  should detect each type of parameter okay
# parameters need to be a vector of names
par.p <-sapply(parameters,output_plots,samp.df = samp.df,rhat=rhat,species_names = species,simple = T)

ggarrange(plotlist=par.p$alpha.p,nrow=4,ncol=1)
ggarrange(plotlist=par.p$dtype1.p,nrow=4,ncol=1)
ggarrange(plotlist=par.p$dtype2.p,nrow=4,ncol=1)
ggarrange(plotlist=par.p$dtype3.p,nrow=4,ncol=1)
ggarrange(plotlist=par.p$mu.beta,nrow=4,ncol=1)
ggarrange(plotlist=sapply(par.p$beta,function(x) x[[6]],simplify = F),nrow=4,ncol=1)

names(par.p$beta) <- covs

ggarrange(plotlist=par.p$alpha.p,nrow=2,ncol=2)
# save plots
wid = 6
hei = 3



# plots of chains, density and parameter estimate 

m.beta <- ggpubr::ggarrange(plotlist = par.p$mu.beta,nrow=2,ncol=1)
ggsave(paste0("effect_plot/Main Effects.1",'_',tribe,group,".png"),m.beta[[1]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.2",'_',tribe,group,".png"),m.beta[[2]],width = wid,height = hei)
ggsave(paste0("effect_plot/Main Effects.3",'_',tribe,group,".png"),m.beta[[3]],width = wid,height = hei)
m.beta

# individual species plots 
ind_plots <- list()

for(i in 1:n.cov){
  ind_plots[[i]] <-  plot_effects(covs[[i]],species=species,cov=i,rhat=rhat,out=out)
  ggsave(paste0( "effect_plot/",covs[[i]],'_',tribe,"_",rq,".png"),ind_plots[[i]],width = wid,height = hei)
}

ind_plots

# check correlation between parameter estimates
cor_plot <- ggs_pairs(samp.df,family="mu.beta" ,lower = list(continuous = "density",alpha=0.2))
ggsave('check_plot/mu_cor_plots.png',cor_plot,width=12,height=9)


#load("jas_out/3_bee_C.3_ID_9.rdata")
#save(out,file="jas_out/3_bee_all_C.3_run.rdata")

#################################################################
############# PP Checks #########################################
#################################################################
# State Model
# aggregate observations where a species was observed at a site


# simulation
covs <- c(base) # covariate list


#  multiple species check
nrep <- 200 #dim(out$sims.list$beta)[1]
beta <- out$sims.list$beta[1:nrep,,] # beta coefficient
gamma <-  out$sims.list$gamma[1:nrep,] # colonization
init <-  out$sims.list$init.occ[1:nrep,] # initial occupancy
alpha <- out$sims.list$alpha.phi[1:nrep,] # intercepts
region.psi <- out$sims.list$region.psi[1:nrep,] # region re
region.cov <- jas_data$nregion

# estimates of population occupancy for species
psi <- sim.pop(psi.start = init,
                gamma = gamma,
                beta = beta,
                alpha = alpha,
                psi =  array(dim=c(nrep,nspecies,sites,time)),
                cov.array = f_array(covs),
                species=T,
               region = region.psi,
               region.cov = jas_data$region,
               complex=T)

# predict occupancy state - write out occupancies to save memory
dir.create("binary_occupancy")
for(i in 1:nrep){
  
  a.x  <-  apply(psi[i,,,],c(2,3),function(x){rbinom(length(x),1,x)})
  saveRDS(a.x,file=paste0("binary_occupancy/",data,"_rep_",i,"_",mod,".rds"))
  
}

# observed occupancy
z <- jas_data[[3]]
z <- z[keep_names,,]

# observation model
out <- occ_output
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
Year <-  observed$Year
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
  
  # could look at transitions states (AHM 2 p243)?
  s_x_e[i] <- sum((colSums(y) - colSums(p))^2/(colSums(p)+e)) # chi square discrepancy over species
  s_x_o[i] <- sum((colSums(occ)  - colSums(p))^2/ (colSums(p)+e))
  v_x_e[i] <- sum((rowSums(y) - rowSums(p))^2/(rowSums(p)+e)) # chi square discrepancy over visits
  v_x_o[i] <- sum((rowSums(occ)  - rowSums(p))^2/ (rowSums(p)+e))
  
  
}


# look at predicted values vs. real values
mean(s_x_e>s_x_o)
mean(v_x_e>v_x_o)

xys  <- seq(min(c(s_x_e,s_x_o))-5, max(c(s_x_e,s_x_o))+5)
xyv  <- seq(min(c(v_x_e,v_x_o))-5, max(c(v_x_e,v_x_o))+5)

mins <- min(c(s_x_e,s_x_o))
ggplot() + geom_point(aes(v_x_e,v_x_o )) +geom_line(aes(xyv,xyv))
ggplot() + geom_point(aes(s_x_e,s_x_o )) +geom_line(aes(xys,xys))

# put lists into date frames
vis_rep. <- do.call(rbind, vis_rep)
sum_rep. <- do.call(rbind, sum_rep)
spe_rep. <- do.call(rbind, spe_rep)

spe_rep.$obs_count <- as.numeric(spe_rep.$obs_count)


act <-colSums( occ)



# plot all species fit
nm <- unique(spe_rep.$snames)
names(nm) <- unique(spe_rep.$snames)
splots <- lapply(nm,function(x){
  
  ggplot() + geom_histogram(data=spe_rep.[spe_rep.$snames==x,] ,aes(x=obs_count))+
    geom_vline(aes(xintercept=act[names(act)==x]),color="blue")+ggtitle(paste(x))
})

indplots <-ggarrange(plotlist=splots,ncol=4,nrow=4)

for(i in 1:length(indplots)){
  
  ggsave(paste0("model_fit/",tribe,"_",mod,"_",data,"_species_fit_",i,".png"), indplots[[i]],width=10,height=8)
  
}

# overall model fit
sum_p <- ggplot() + geom_histogram(data= sum_rep.
                                   ,aes(x=y_sum),bins =40)+
  geom_vline(aes(xintercept=sum(rowSums(occ))),color="blue")+ggtitle(paste0('Total Count Estimates \nModel = ',data," ",mod))

ggsave(paste0("model_fit/",tribe,"_",mod,"_",data,"_all_fit.png"), sum_p,width=4,height=4 )


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


os <-readRDS("osmia.rds")
sam <- ggs(os$BUGSoutput)
unique(sam$Parameter)
plot(os$BUGSoutput)
traceplot(os$BUGSoutput$sims.list$dtype2.p)
s<-os$BUGSoutput$programe
