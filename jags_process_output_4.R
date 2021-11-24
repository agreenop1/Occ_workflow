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

# iter.index = UID in file path needs to be related to unique index containing different chunks of the chain (iteration numbers)
# chain.index = if comb.chain = T chains need to be combined. CID used to identify the different chain in the file path
# file.path = e.g. paste0("Model_outputs/Mod_sumRQ_500_C.","CID","_ID_","UID",".rdata")
# by.it = number of increases in iterations in each file
# iterations = iterations to be used
# verbose = T/F output to console



source("Occ_workflow/combine_chain_4.1.R")

file.path = paste0("Jasmin_outputs/", "spi.B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("beta","deviance","gamma",
                  "init.occ"),iter.index=1,chain.index=2:3,summary=T,file.path=file.path,by.it=1000,iterations=1000,verbose=T)

# output <- summary_chains(out,comb.chain = T,keep.samples = F)

# model 
print(out)
samples = out$samples
samp.df <- ggmcmc::ggs(samples)

# co-variates in order used in model
covs <- c("Temperature","RQsum","Semi.natural")

# species in order they are in occ dataframe!
jas_data <- readRDS("Model_data/data_spiders_all.499_1994.2010.rds")
nspecies <- colnames(jas_data[[1]][-1])

# parameters of interest
params_jname <- list("beta","mu.beta","deviance","gamma","init.occ")

par <- paste0(params_jname,collapse = "|")
vars. <- varnames(out$samples)
vars <- vars.[stringr::str_detect(vars.,par)]
#
plots <- function(v1){
tr <- ggs_traceplot(samp.df[samp.df$Parameter==v1,])+
  ggtitle(paste(cnames[c],"-",rnames[r],"- rhat =",round(out$Rhat$beta[r,c],2)))
ds <- ggs_density(samp.df[samp.df$Parameter==v1,])
coef <- ggs_caterpillar(samp.df[samp.df$Parameter==v1,])+
  geom_vline(xintercept = 0,linetype="dashed")

ggpubr:: ggarrange(tr,ds,coef)
}
# convergence plots & rhat statistics based on selected parameters #
####################################################################
# parameters indexed in matrix e.g. coefficient x species [1,2]
# rownames in order
rnames = covs; row.it = length(covs)
cnames = species; col.it = length(species)

par.plots <- list()

for (c in 1:col.it){
  
  cov <- list()
  
  for (r in 1:row.it){
  v1 <- gsub("r",r,paste0("beta[r,c]"))
  v1 <- gsub("c",c,v1)
  
  cov[[r]] <- plots(v1)
  names(cov)[r] <- rnames[r]
  }

 gm <- paste0("gamma[",c,"]")    
 io <- paste0("init.occ[",c,"]") 
  cov2 <- list(gamma= plots(gm),psi=plots(io))
par.plots[[c]] <- c(cov,cov2)
names(par.plots)[c] <- cnames[c]
}

################################################## 
hpara <- c("mu.beta","tau.beta")
var.plots <- list()

for (r in 1:row.it){
  
  cov <- list()
  
  for (i in 1:length(hpara)){
    v <- hpara[i]
    v1 <- paste0(v,"[",r,"]")
    
    tr <- ggs_traceplot(samp.df[samp.df$Parameter==v1,])+
      ggtitle(paste(rnames[r],"- rhat =",round(out$Rhat[[v]][r],2)))
    ds <- ggs_density(samp.df[samp.df$Parameter==v1,])
    coef <- ggs_caterpillar(samp.df[samp.df$Parameter==v1,])+
      geom_vline(xintercept = 0,linetype="dashed")
    
    p1 <- ggpubr:: ggarrange(tr,ds,coef)
    cov[[i]] <- p1
    names(cov)[i] <- hpara[i]
  }
  
 var.plots[[r]] <- cov
 names(var.plots)[r] <- rnames[r] 
}
var.plots
##
x <- names(out$Rhat)

rhat <- function(x){
cat(x,"\n")
x1 <- out$Rhat[[x]]
if(is.matrix(x1)){

  x1[x1<1.05] <- "Converged (Rhat < 1.05)"
  x1[x1!="Converged (Rhat < 1.05)"] <- "Not converged (Rhat > 1.05)"
  x1 <- t(x1)
  rownames(x1) <- species
  colnames(x1) <- covs
  x1
}else{
  x1[x1<1.05] <- "Converged (Rhat < 1.05)"
  x1[x1!="Converged (Rhat < 1.05)"] <- "Not converged (Rhat > 1.05)"
  
  if (length(x1)==length(species)){ names(x1) <- species}
  if (length(x1)==length(covs)){ names(x1) <- covs}
  x1
 }
}
c<- lapply(x,rhat)
names(c) <- x
c
#################################################################
############# Mean coefficient for each species #################
#################################################################


plot_effects <- function(x,species){
  pdat <- data.frame(matrix (nrow=length(species),ncol=4))
  colnames(pdat) <- c("Species","Mean","Lb","Ub")
  pdat[,1] <- species
  pdat[,2] <- means
  pdat[,3] <- lb
  pdat[,4] <- ub
  
  
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

pr <- 1
means <- out$mean$beta[pr,]
lb <- out$q2.5$beta[pr,]
ub <- out$q97.5$beta[pr,]
plot_effects("temp",species=species)

