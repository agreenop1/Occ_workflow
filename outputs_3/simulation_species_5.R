################################################################################
################################################################################
##################### simulate different scenarios #############################
library(brms)

source("Occ_workflow_V2/outputs_3/combine_chain_4.1.R")

# read in parameters
para <- c("mu.beta","gamma","beta","init.occ","alpha.phi","mu.gamma")

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
iterations. <- 300

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2016.rds"))


para = "phi"
# compile the chains from the model
out <- comb_daisy(parameters=para,
                  iter.index=18:18,chain.index=1:1,summary=F,file.path=file.path,by.it=500,
                  it.used=its,
                  iterations=1000,verbose=T)
simulation.l(out$samples)

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

plots_trend <- list()
slope_plots <- list()
################################################################################
############################### Simulation #####################################
################################################################################
for(s in 1:94){

names  <- colnames(jas_data[[1]][-1])[s]

# out latent occupancy
beta <- out$sims.list$beta[,,s]

# gamma
gamma <-  out$sims.list$gamma[,s]
summary(mean(gamma)) 

# initial occupancy
init <- out$sims.list$init.occ[,s]
summary(mean(init)) 

# alpha.phi - intercept
alpha.phi <- out$sims.list$alpha.phi[,s]
summary(mean(alpha.phi)) 
  
# covariates
covars = jas_data[[4]]
base <- list(temp.m=covars$mean_temp,temp.a=covars$temp_anom,semi=covars$semi,agri=covars$agri,RQA=covars$RQsum_A,RQM=covars$RQsum_M)
covs <- c(base) 
covs1 <- c(base)



# min vlue
min(covs$RQA)
covs1$RQA[,1:12] <- -22.59073

min(covs$RQM)
covs1$RQM[,1:12] <- -5.855139


###############################################################
##################### Run simulations #########################
###############################################################
# parameters for simulation
psi1 <- sim.pop(psi.start = init,
                gamma = gamma,
                alpha = alpha.phi,
                psi = array(dim = c(iterations.,sites,time)),
                cov.array = f_array(covs))

psi2 <- sim.pop(psi.start = init,
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
plots_trend[[names]] <- ggplot(data=m)+geom_line(aes(x=year,y=mean,color=name))+
         geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=name),alpha=0.2)+ylim(0,1)+
         ggtitle(names)
  

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

diff <- s1-s2

m2 <- summ(diff,year = seq(1996,2018,2)) 

slope_plots[[names]] <- ggplot(data=m2)+geom_point(aes(x=year,y=mean))+
  geom_errorbar(aes(x=year,ymin=lb,ymax=ub),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  ggtitle(names)
}
