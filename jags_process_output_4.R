library(BRCmap)
library(dplyr)
library(ggplot2)
library(jagsUI)
library(forcats)
library(dplyr)

# iter.index = UID in file path needs to be related to unique index containing different chunks of the chain (iteration numbers)
# chain.index = if comb.chain = T chains need to be combined. CID used to identify the different chain in the file path
# file.path = e.g. paste0("Model_outputs/Mod_sumRQ_500_C.","CID","_ID_","UID",".rdata")
# by.it = number of increases in iterations in each file
# thin = thin rate
# iterations = iterations to be used
# comb.chain = T/F combining chains run in parallel
# verbose = T/F output to console
library(jagsUI)

source("Process_output/combine_chain.R")

file.path = paste0("Jasmin_outputs/", "bee.B.all_C.","CID","_run.rdata")
out <- comb_daisy(parameters=c("beta","deviance"),iter.index=1,chain.index=1:3,file.path=file.path,by.it=500,iterations=500,verbose=T)
output <- summary_chains(out,comb.chain = T,keep.samples = T)
plot(output,"mu.beta[1]")

print(output)


# mod output
mod <- readRDS("ex.mod.rds")

# species level parameters
eco_params <- list("beta1","beta2")

# site x year parameters
obs_params <- list("z","phi","psi")

species <- c("a","b","c") # in order they are in occ dataframe!


####################################################################
# convergence plots & rhat statistics based on selected parameters #
####################################################################
eco_conv_plot <- lapply(eco_params,function(x,y){ recordPlot(plot(y,x))}, y=mod)
obs_conv_plot <- lapply(obs_params,function(x,y){ recordPlot(plot(y,x))}, y=mod)

# check convergence
rhat <- mod$Rhat

rhat_asses <- function(x,rhat,species, M=F){
  para <- x[[1]]
  if (!M){y <- data.frame(matrix(ncol=length(x),nrow=length(species)))
  } else {y <-  array(dim=c(length(x),ncol(rhat[[para]][1,,])+1,length(species)))}
  
  if(!M){
    for (i in 1:ncol(y)){
      y[1:nrow(y),i]  <- round(rhat[x[[i]]][[1]],3) 
      colnames(y)[i] <- x[[i]]
    }
    
    y[,ncol(y)+1]<- apply(y,1,function(x){
      if(!is.numeric(x)) stop("Not all parameters are numeric")
      if(all(x<1.05)) "All params converged (rhat < 1.05)"else "At least 1 param did NOT converge (rhat >1.05)"})
    
    colnames(y)[ncol(y)] <- "Convergence status"
    rownames(y)[1:length(species)] <- species
    y
    
  }else{
    for (i in 1:length(species)){
      for (k in 1:length(x)){
        
        vals <- rhat[[x[[k]]]][i,,]
        
        y[k,,i] <- c(x[[k]],round(apply(vals,2,function(x)sum(x<1.05,na.rm=T))/nrow(vals),4)*100)
        
      }
    }
    y
  }
}


# ecology & observation mod rhat statistics

rhat_eco <- rhat_asses(eco_params,rhat,species,M=F)
rhat_s <- rhat_asses(obs_params,rhat,species,M=T)

#################################################################
############# Mean coefficient for each species #################
#################################################################
means <- mod$mean
lb <- mod$q2.5
ub <- mod$q97.5

plot_effects <- function(x,species){
  pdat <- data.frame(matrix (nrow=length(species),ncol=4))
  colnames(pdat) <- c("Species","Mean","Lb","Ub")
  pdat[,1] <- species
  pdat[,2] <- means[[x]]
  pdat[,3] <- lb[[x]]
  pdat[,4] <- ub[[x]]
  
  
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

lapply(eco_params,plot_effects,species=species)


##########################################################################
############### Parameter effects on species persistence #################
##########################################################################
sims <- mod$sims.list

# covariate data
chems <- read.csv("Weight_applied_byChemgroup.csv")
chems <- split(chems,chems$chemgroup)

x<- chems[[1]]

# parameter values
FUNG <- seq(min(chems[[1]]$weight_applied),max(chems[[1]]$weight_applied),length.out=100)-mean(chems[[1]]$weight_applied)
FUIN <- seq(min(chems[[2]]$weight_applied),max(chems[[2]]$weight_applied),length.out=100)-mean(chems[[2]]$weight_applied)
HERB <- seq(min(chems[[3]]$weight_applied),max(chems[[3]]$weight_applied),length.out=100)-mean(chems[[3]]$weight_applied)
INSE <- seq(min(chems[[4]]$weight_applied),max(chems[[4]]$weight_applied),length.out=100)-mean(chems[[4]]$weight_applied)
MOLL <- seq(min(chems[[5]]$weight_applied),max(chems[[5]]$weight_applied),length.out=100)-mean(chems[[5]]$weight_applied)


## values for forming predictions
vari <- data.frame(FUIN, INSE)
mean_v <- c(mean(chems[[2]]$weight_applied),mean(chems[[4]]$weight_applied))
betas <- list("beta1","beta2")
intercept <- "alpha.phi"

spec_vals <- list()
para_vals <- list()

for (i in 1:ncol(vari)){
  for (k in 1:length(species)){
    beta <- betas[[i]]
    vals <- sapply(vari[,i],function(x) plogis(sims[[intercept]][,k]+ sims[[beta]][,k]*x))
    means <- colMeans(vals)
    lb <- apply(vals,2,quantile,probs = 0.025)
    ub <- apply(vals,2,quantile,probs = 0.975)
    spec_vals[[k]] <- data.frame(means,lb,ub,vari[,i]+mean_v[[i]],colnames(vari[i]),species[[k]])
  }
  para_vals[[i]] <- spec_vals
}


spec_effs <- lapply(para_vals,
                    function(x)lapply(x,
                                      function(x) ggplot(data=data.frame(x)) +  
                                        geom_ribbon(aes(ymin=ub,ymax=lb,x=x[,4]),alpha=0.4)+ 
                                        geom_line(aes(x=x[,4],y=means))+
                                        xlab(paste(x[1,5],"weight"))+ylab("Persistence")))

spec_effs[[2]]

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

