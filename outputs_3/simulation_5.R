################################################################################
################################################################################
##################### simulate different scenarios #############################
library(brms)
library(BRCmap)
library(tidyverse)
library(dplyr)
source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")
UK <-  readRDS("UK_map.rds")
# read in parameters
para <- c("mu.beta","gamma","beta","init.occ","mu.alpha.phi","mu.gamma")

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
iterations. <- 1800

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))

# compile the chains from the model
#out <- comb_daisy(parameters=para,
#                  iter.index=18:18,chain.index=1:3,summary=T,file.path=file.path,by.it=500,
#                  it.used=its,
#                  iterations=1000,verbose=T)

out <- readRDS("jas_out/summary_p.rds")[[2]]


# out latent occupancy
beta <- out$sims.list$mu.beta
colMeans(beta)
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
  
# covariates
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)



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
sim.pop <-  function(psi,psi.start,alpha,gamma,cov.array){

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


# summarize output
summ <- function(x,year){
  m1 <- apply(x,c(3),mean)
  m2 <- apply(x,c(3),quantile,probs = 0.025)
  m3 <- apply(x,c(3),quantile,probs = 0.975)
  
  trend <- data.frame(mean=m1,lb=m2,ub=m3)
  trend$year <- year 
  trend}
###############################################################
##################### Run simulations #########################
###############################################################
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
