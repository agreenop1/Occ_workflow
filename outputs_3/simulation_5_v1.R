################################################################################
################################################################################
##################### simulate different scenarios #############################
# look at spacial extent of declines and the severity
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
  m1 <- apply(x,c(3),mean,na.rm=T)
  m2 <- apply(x,c(3),quantile,probs = 0.025,na.rm=T)
  m3 <- apply(x,c(3),quantile,probs = 0.975,na.rm=T)
  m4 <-round( apply(x,c(2,3),quantile,probs = 0.05,na.rm=T),3)
  m5 <-round( apply(x,c(2,3),quantile,probs = 0.95,na.rm=T),3)
  m6 <-round( apply(x,c(2,3),sd,na.rm=T),3)
  
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


sp1<- ggplot(data=m3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")

ggsave("effect_plot/sp1_time.png",plot= sp1,width = 5,height=3)

###############################################################
# Change Maps 
# summarize output
summs <- function(x,year){
  m1 <-round( apply(x,c(2,3),mean,na.rm=T),3)
  m2 <-round( apply(x,c(2,3),quantile,probs = 0.025,na.rm=T),3)
  m3 <-round( apply(x,c(2,3),quantile,probs = 0.975,na.rm=T),3)
  m4 <-round( apply(x,c(2,3),quantile,probs = 0.05,na.rm=T),3)
  m5 <-round( apply(x,c(2,3),quantile,probs = 0.95,na.rm=T),3)
  m6 <-round( apply(x,c(2,3),sd,na.rm=T),3)
if(is.null(year)){
list(mean = m1,lb = m2, ub = m3,lb5=m4,ub90=m5)}else{
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
                               aes(x = E, y = N, fill =sd))+
  scale_fill_continuous(type = "viridis", name = "")+ theme( axis.text.y = element_text(size=10),
                                                             axis.text.x=  element_text(size=10),
                                                             axis.title.x= element_text(size=11), 
                                                             axis.title.y =element_text(size=11),panel.grid.major = element_blank(), 
                                                                                    panel.grid.minor = element_blank(),
                                                                                    panel.background = element_blank(),
                                                           )+
  ggtitle("SD change in occupancy")
)
plots_map[[1]]
ggsave("spatial_plot.png",dpi=2000,plot= plots_map[[1]],width = 4.5,height=4.5)



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

m1 <- summ(psi1,year=seq(1994,2018,2))
m1$name <- "With pesticide"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Without pesticide"

m <- rbind(m1,m2)


#diffp<-  psi1-psi2
#n_site <- apply(diffp,c(1,2,3),function(x){if_else(x==0,x,1)})
#n_site_s <- apply(n_site,c(1,3),sum) 
##difference <-  apply(diffp,c(1,2,3),function(x) ifelse(x==0,NA,x) )

diffp <-  psi1-psi2
slope <- function(x){
  l <- dim(x)[1]
  time <- dim(x)[3]
  
  out <- array(dim=c(l,dim(x)[2],time-1))
  
  for(i in 1:l){
    
    out[i,,]  <- (x[i,,-1]-x[i,,1:12])/2
    
  }
  out
}

sdif <- slope(diffp)
m3<- summ(sdif,year =seq(1996,2018,2))
library(ggplot2)
ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name),size=1)+
  geom_ribbon(aes(x=year,ymin=lb5,ymax=ub90,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")


m3 <- summ(diffp,year = seq(1994,2018,2)) 

temp_ef<- ggplot(data=m3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")
ggsave("effect_plot/temp_ef.png",plot= temp_ef,width = 5,height=3)
  ###############################################################
# site level changes between the first and the last year
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


tmaps<- ggarrange(plotlist=plots,nrow=2,ncol=3)
map<- ggarrange(plotlist=list(plots[[1]],plots[[11]],plots[[12]]) ,nrow=1,ncol=3)
ggsave("temp_plot1.png",plot=  tmaps$`1`,width = 10,height=7)
ggsave("temp_plot2.png",plot= tmaps$`2`,width = 10,height=7)
ggsave("temporal_rq.png",plot= map,width = 16,height=7,dpi = 2000)
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

############################################################################
##################### Persistence (%) Marginal Effects #########################
############################################################################
covsr <- readRDS("covars_wide.rds")

# temperature marginal effect ##############################################
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)
temps <- covsr$temp_anom[-1]+covsr$mean_temp[-1]
covs <- c(base) 
covs1 <- c(base)
me <-mean(as.matrix( covsr$mean_temp[-1]))
mn<- min(covsr$mean_temp[-1])-me
mx<- max(covsr$mean_temp[-1])-me
alpha = alpha.phi
c1 <- seq(mn,mx, length.out =100)
phi <-  apply( t(matrix(ncol=iterations.,nrow=100,(c1))) *beta[,1] ,2,function(x){inv_logit_scaled( x+alpha)})
summp <- function(x){
  m1 <- apply(x,c(2),mean,na.rm=T)
  m2 <- apply(x,c(2),quantile,probs = 0.025,na.rm=T)
  m3 <- apply(x,c(2),quantile,probs = 0.975,na.rm=T)

  
  trend <- data.frame(mean=m1,lb=m2,ub=m3)

  trend}
b1 <- summp(phi)

tmp <- ggplot(data=b1)+geom_line(aes(x=c1+me,y=mean*100),size=1)+
  geom_ribbon(aes(x=c1+me,ymin=lb*100,ymax=ub*100),alpha=0.2)+
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

# risk normally marginal effect ##############################################
rq <- covsr$RQsum_A[-1]+covsr$RQsum_M[-1]
mx=max(rq)
mn=min(rq)-mean(as.matrix(rq))
me=mean(as.matrix(rq))
sdrq=sd(as.matrix(rq))
quantile(as.matrix(covsr$RQsum_A[-1]),1)
quantile(as.matrix(covsr$RQsum_A[-1]),0.95)
quantile(((rq/covsr$RQsum_M[-1])/covsr$RQsum_M[-1])*100,0.95,na.rm=T)
max(as.matrix(covsr$RQsum_A[-1]))
phi <-  apply( t(matrix(ncol=iterations.,nrow=100,(seq(mn,mx, length.out =100)))) *beta[,5] ,2,function(x){inv_logit_scaled( x+alpha)})

b1 <- summp(phi)



ggplot(data=b1)+geom_line(aes(x=seq(mn,mx, length.out =100),y=mean),size=1)+
  geom_ribbon(aes(x=seq(mn,mx, length.out =100),ymin=lb,ymax=ub),alpha=0.2)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Persistence (%)")+xlab("Risk Anomaly")


# risk normally marginal effect ##############################################
rq <- covsr$RQsum_A[-1]+covsr$RQsum_M[-1]
((rq/covsr$RQsum_M[-1])/covsr$RQsum_M[-1])*100
me=mean(as.matrix(rq))
mx=60.01417
mn=min(rq)-mean(as.matrix(rq))
max=max(rq)
min=min(rq)

phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,5] ,2,function(x){inv_logit_scaled( x+alpha)})

per<-(((seq(mn,mx, length.out =200)+me)-me)/me)*100

b1 <- summp(phi)
nf<- quantile(((rq/covsr$RQsum_M[-1])/covsr$RQsum_M[-1])*100,0.95,na.rm=T)
b2 <- summp(data.frame( inv_logit_scaled(alpha+( 42.71109*beta[,5]))))
(sdrq+me)-me
b3 <- summp(data.frame( inv_logit_scaled(alpha+(60.01417*beta[,5]))))
b4 <- summp(data.frame( inv_logit_scaled(alpha)))



me=mean(as.matrix(rq))
covsr$RQsum_A %>%filter_all(any_vars(. %in% c(max(covsr$RQsum_A[-1]))))

rqap<-ggplot(data=b1)+geom_line(aes(x=per,y=mean*100),size=1)+
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
rq <- covsr$semi[-1]
mx=max(rq)-mean(as.matrix(rq))
mn=min(rq)-mean(as.matrix(rq))
max=max(rq)
min=min(rq)
phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,3] ,2,function(x){inv_logit_scaled( x+alpha)})

b1 <- summp(phi)
b2 <- summp(data.frame( inv_logit_scaled(alpha)))


me=mean(as.matrix(rq))





semi <- ggplot(data=b1)+geom_line(aes(x=seq(min,max, length.out =200),y=mean*100),size=1)+
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
ggsave("semi_eff.png",plot= semi,width = 6,height=3)


# agriculture marginal effect #######################################

rq <- covsr$semi[-1]
mx=max(rq)-mean(as.matrix(rq))
mn=min(rq)-mean(as.matrix(rq))
max=max(rq)
min=min(rq)
phi <-  apply( t(matrix(ncol=iterations.,nrow=200,(seq(mn,mx, length.out =200)))) *beta[,3] ,2,function(x){inv_logit_scaled( x+alpha)})
b1 <- summp(phi)
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

################################################################################
# PP check #####################################################################
total_obs <- cbind(occ=rowSums( jas_data[[1]][-1]),jas_data[[2]])%>%
               group_by(site_5km,TP)%>%summarise(total=sum(occ))
plot(density(log(total_obs$total+0.1)))
