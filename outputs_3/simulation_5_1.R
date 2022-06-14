################################################################################
################################################################################
##################### simulate different scenarios #############################
library(brms)

source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")

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
iterations. <- 1500

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))


out <- readRDS("jas_out/summary_p.rds")[[2]]


# out latent occupancy
beta <- out$sims.list$mu.beta
colMeans(beta)
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
covs <- c(base) 
covs1 <- c(base)



# min vlue
min(covs$RQA)
covs1$RQA[,1:12] <- -22.59073 # e.g 0 pesticide applied

min(covs$RQM)
covs1$RQM[,1:12] <- -5.855139


# covariates
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

library(ggplot2)
tp  <- ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name),size=1)+
         geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
  )+ylab("Occupancy")+xlab("Year")

tp$labels$fill <- "Scenario"
tp$labels$colour <- "Scenario"
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
diff<-  psi1-psi2
m2 <- summ(diff_s,year = seq(1996,2018,2)) 
m3 <- summ(diff,year = seq(1994,2018,2)) 

s<-ggplot(data=m2)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in slopes")+xlab("Year")

o  <-ggplot(data=m3)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Difference in occupancy")+xlab("Year")

ggpubr::ggarrange(tp,s,o)
###############################################################################
######################## Occupancy check ######################################
###############################################################################
s=1
gamma <-     out[[1]]$sims.list$gamma[,s]
alpha.phi <- out[[1]]$sims.list$alpha.phi[,s]
init <-   out[[1]]$sims.list$init.occ[,s]
beta <- out[[1]]$sims.list$beta[,,s]

psi2 <- sim.pop(psi.start = init,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs))


# calculate psi

psi = array(dim = dim(phi) )

psi[,1:sites,1] <- init

phi <- out[[2]]$phi[,s,,]

for(i in 1:sites){
  for(t in 2:time){
  
    psi[,i,t] <-  (psi[,i,t-1]* phi[,i,t]) + (1- psi[,i,t-1])*gamma
  }
}

all(psi==psi2)

m1 <- summ(psi,year=seq(1994,2018,2))
m1$name <- "Actual"
m2 <- summ(psi2,year = seq(1994,2018,2)) 
m2$name <- "Predicted"

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

