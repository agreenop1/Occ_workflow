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

file="lad" 
data="ladybirds"

source("Occ_workflow/combine_chain_4.1.R")

file.path = paste0("Jasmin_outputs/",file ,".B.all_C.CID_ID_UID.rdata")
out <- comb_daisy(parameters=c("alpha.phi","mu.alpha.phi","tau.alpha.phi","alpha.p",
                               "dtype1.p","dtype2.p","dtype3.p",
                               "mu.d3.p","tau.p3","mu.d1.p","tau.p1","mu.d2.p","tau.p2","tau.p",
                               "gamma","mu.gamma","tau.gamma",
                               "beta","mu.beta","tau.beta",
                               "init.occ"),
                  iter.index=1:32,chain.index=1:3,summary=T,file.path=file.path,by.it=1000,iterations=10000,verbose=T)

# output <- summary_chains(out,comb.chain = T,keep.samples = F)

# model 
print(out)
samples = out$samples
samp.df <- ggmcmc::ggs(samples)

# co-variates in order used in model
covs <- c("Temperature_anom","Temp_mean","RQsum","Semi.natural")

# species in order they are in occ dataframe!
jas_data <- readRDS(paste0("Model_data/data_",data,"_all.499_1994.2010.rds"))
species <- colnames(jas_data[[1]][-1])

nspecies <- length(species)
# parameters of interest
params_jname <- list("beta","mu.beta","deviance","gamma","init.occ")

par <- paste0(params_jname,collapse = "|")
vars. <- varnames(out$samples)
vars <- vars.[stringr::str_detect(vars.,par)]
#
plots <- function(v1){
tr <- ggs_traceplot(samp.df[samp.df$Parameter==v1,])
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
 al <- paste0("alpha.phi[",c,"]") 
 p <- paste0("dtype1.p[",c,"]") 
  cov2 <- list(gamma= plots(gm),psi=plots(io),a.phi=plots(al),p=plots(p))
par.plots[[c]] <- c(cov,cov2)
names(par.plots)[c] <- cnames[c]
}
par.plots[[9]]$p
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

pr <- 3
means <- out$mean$beta[pr,]
lb <- out$q2.5$beta[pr,]
ub <- out$q97.5$beta[pr,]
plot_effects("temp",species=species)
par.plots$`Subcoccinella vigintiquattuorpunctata`
save.image("dtype1.1.jpeg")
jpeg("dtype1_4.JPG")
save

for (i in 1:length(species)){print(par.plots[[i]]$a.phi)}
plot(out,"alpha.phi")

