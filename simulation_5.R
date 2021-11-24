#########################################################################
# 5. Code for simulation of dynamics bayesian occupancy models          #
#########################################################################
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
plot_div<-function(x,title,subtitle,x.axis.texts,x.title.texts,color,y.axis.texts){
  
  library(ggplot2)
  
  # figure for diversity measures
  # standardised effect size plot
  # trend plot
  fig<-ggplot(x)+geom_line(aes(y=Mean,x=Year),color=color)+
    geom_errorbar(aes(x=Year,ymin=LB_95, ymax=UB_95,width=0.6),color=color)+
    theme(axis.text.y = element_text(size=7),
          axis.text.x= x.axis.texts ,
          axis.title.x=x.title.texts, 
          axis.title.y = y.axis.texts,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title=element_text(size=8), 
          legend.text=element_text(size=7),
          legend.key.size = unit(0.5,"line"),
          legend.position="none",
          plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(face="bold",size=9),
          plot.subtitle=element_text(size=9))+ylab("Occupancy")+ggtitle(title,subtitle = subtitle)+
    geom_segment(aes(x=Year,xend=Year,y=LB_80, yend=UB_80),color=color, size = 2)+
    geom_point(aes(y=Mean,x=Year), size = 1)
  
  # list of plots    
  fig
}


source("Occ_workflow/combine_chain_4.1.R")

file.path = paste0("Jasmin_outputs/", "lad.B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("beta","deviance","gamma",
                               "init.occ"),iter.index=1:40,chain.index=1:3,summary=T,file.path=file.path,by.it=500,iterations=10000,verbose=T)

# species in order they are in occ dataframe!
jas_data <- readRDS("Model_data/data_ladybirds_all.499_1994.2010.rds")

# output <- summary_chains(out,comb.chain = T,keep.samples = F)

# model 
print(out)



#######################################################################################
# data in format from jags data prep
occup<- jas_data[[1]] #occ data
visit<- jas_data[[2]]# visit info
zobs<- jas_data[[3]]# init values
closure.period <- jas_data[[5]] #closure period

# check sample sizes in original data set
nsites <- length(unique(visit$site_5km)) # sites
nclosure <- length(unique(visit$TP)) # number of time points
nspecies <- ncol(occup[-1]) # num. species
nits <- out$mcmc.info$n.samples

out$model[[1]][[1]]$model

# covariate formula
covars <- jas_data[[4]]# covariates
covs <- list(temp=covars$temp_anom,rqsum=covars$RQsum,semi=covars$semi) # named list
inters <- NULL # can supply a vector of which variables have interactions based on position in above list
###################################################################
# covarite processing

# create array for covs
row = nrow(covars[[1]]); col = ncol(covars[[1]]); nmat = length(covs)
cov.array <- array(dim=c(row,col,nmat))

for(n in 1:nmat){
  cov.array[,,n]  <- as.matrix(covs[[n]])
}

# formula for main effects
formula <- paste0("beta[",1:nmat,",s]*COVS[i,t-1,",1:nmat,"]",collapse ="+")

# if any interactions add to formula

if(!is.null(inters)){
  coms <- combn(inters,2)  
  
  for (n in 1:ncol(coms)){
    ints1 <- paste0("+beta[",n+nmat,",s]*COVS[i,t-1,",coms[1,n],"]*COVS[i,t-1,",coms[2,n],"]")
    formula <- paste0(formula,ints1) }
  
  nmat <- ncol(coms)+nmat
}

vars <- strsplit(formula,split="+",fixed=T)[[1]]

# output formula
cat("Ecological model covariate structure","\n")
cat(vars,sep="+ \n")
cat(formula,"\n")

sims <- out$sims.list
out <-NULL
nits = nrow(sims$gamma)
sit=300
psi <- array(dim=c(sit,nspecies,nsites,nclosure))

for (i in 1:nsites){
  psi[,,i,1] <- sims$init.occ[sample(1:nits,sit),]
}
phi<- array(dim=c(sit,nspecies,nsites,nclosure))
gamma <- sims$gamma
alpha.phi <- sims$gamma
beta <- sims$beta
COVS = cov.array


for (s in 1:nspecies){
  for (i in 1:nsites){
    for (t in 2:nclosure){
      
      phi[,s,i,t] <-   plogis (  sample(alpha.phi[,s],sit) +
                                   sample(beta[,1,s],sit)*COVS[i,t-1,1] +
                                   sample(beta[,2,s],sit)*COVS[i,t-1,2] +
                                   sample(beta[,3,s],sit)*COVS[i,t-1,3]) 
      
      psi[,s,i,t] <- psi[,s,i,t-1]*phi[,s,i,t] + (1-psi[,s,i,t-1])*sample(gamma[,s],sit)
      
    }
  } 
}

# 
ex.closure=5; ex.COV = NULL; 
# mean trend
x=psi
s.mean<- apply(x,c(2,4),mean)   
s.025<-  apply(x,c(2,4),quantile,probs = 0.025)
s.0975<- apply(x,c(2,4),quantile,probs = 0.975)
s.010<-  apply(x,c(2,4),quantile,probs = 0.10)
s.090<-  apply(x,c(2,4),quantile,probs = 0.90)
a=list()
# 
for (i in 1:nspecies){
x1=data.frame(Mean=s.mean[i,],
           LB_95  =    s.025[i,], 
           UB_95 =   s.0975[i,],
           LB_80 =    s.010[i,], 
           UB_80 =    s.090[i,],
           Year = seq(1994,2012,2))

a[[i]]=plot_div(x1,title=NULL,subtitle = "Trend", color="#3399FF",x.axis.texts = element_text(size=8),x.title.texts =element_text(size=8),
         y.axis.texts = element_text(size=8))
}
a

# maps
# percent or 0 beta
func.det <- array(dim=c(nclosure,nspecies,nsites,sit))
for(ni in 1:sit ){
  for (s in 1:nspecies){
    for (i in 1:nsites){
      for (t in 1:nclosure){
        
        func.det[t,s,i,ni] <- psi[ni,s,i,t]
        
      }
    }
  }
}


betapart::beta.pair.abund(func.det[,,1,1])

##########################################################################
############### Parameter effects across England #########################
##########################################################################
covar <- readRDS("chems_W_SID.csv")
psi <- mod$sims.list["psi"]

var.x <- "Insecticide"

dat.x  <- covar[[var.x]]
var_mean <- data.frame((colMeans(covar[[var.x]][,-1],na.rm = T)))
colnames(var_mean)[1] <- c("mean")
var_mean$year  <- rownames(var_mean)
var_mean <- var_mean[order(var_mean$mean),]
min.x <- var_mean[1,2]
max.x <- var_mean[nrow(var_mean),2]

x.map <- cbind(dat.x$gr, dat.x[[min.x]]/dat.x[[max.x]])
y.map <- 
  
  
  ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = filter(p, year == "1994"), 
                               aes(x = E, y = N, fill = weight_applied))+
  scale_fill_continuous(type = "viridis", name = "Crop Area (ha)")

proj4string(UK$britain)
c<-UK$britain
gelman.diag(p,multivariate = F)
dic(mod$model)

##############################################
trend<-function(x){
  out<- data.frame (matrix (ncol=6,nrow = 10))
  colnames(out)<-c("Year","Mean","95_LB","95_UB","80_LB","80_UB")
  out[1:10,1]<- seq(1970,2015,by=5)
  out[1:10,2]<- apply(x,1,mean)   
  out[1:10,3]<- apply(x,1,quantile,probs = 0.025)
  out[1:10,4]<- apply(x,1,quantile,probs = 0.975)
  out[1:10,5]<- apply(x,1,quantile,probs = 0.10)
  out[1:10,6]<- apply(x,1,quantile,probs = 0.90)
  out}


 
