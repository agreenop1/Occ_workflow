library(BRCmap)
library(dplyr)
library(ggplot2)
library(jagsUI)
library(forcats)
library(dplyr)
library(ggmcmc)
library(reshape2)
# species in order they are in occ dataframe!
jas_data <- readRDS("Model_data/data_ladybirds_all.499_1994.2010.rds")

##########################################################################################
plot_div<-function(x,title,subtitle,x.axis.texts,x.title.texts,color,y.axis.texts,group=F){
  
  library(ggplot2)
  
  if(group!=T){
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
    fig} else{
      # figure for diversity measures
      # standardised effect size plot
      # trend plot
      dodge <- position_dodge(width=0.5)  
      
      fig<-ggplot(x)+geom_line(aes(y=Mean,x=Year,colour=Cond),position = dodge)+
        geom_errorbar(aes(x=Year,ymin=LB_95, ymax=UB_95,width=1,colour=Cond),position = dodge)+
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
              legend.position="left",
              plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"),
              axis.line = element_line(colour = "black"),
              plot.title = element_text(face="bold",size=9),
              plot.subtitle=element_text(size=9))+ylab("Occupancy")+ggtitle(title,subtitle = subtitle)+
        geom_linerange(aes(x=Year,ymin=LB_80, ymax=UB_80,colour=Cond),size=2,position = dodge)+
        geom_point(aes(y=Mean,x=Year,colour=Cond),position = dodge, size = 1)
      
      # list of plots    
      fig
    }
}


#######################################################################################
# data in format from jags data prep
occup<- jas_data[[1]] #occ data
visit<- jas_data[[2]]# visit info


# check sample sizes in original data set
nsites <- length(unique(visit$site_5km)) # sites
nclosure <- length(unique(visit$TP)) # number of time points
nspecies <- ncol(occup[-1]) # num. species
species <- colnames(occup[-1])
site_id <- jas_data[[8]]
rm(jas_data)
# mean trend


sumr <- function(x,y=c(2,4)){list(
  s.mean<- apply(x,y,mean)   ,
  s.025<-  apply(x,y,quantile,probs = 0.025),
  s.0975<- apply(x,y,quantile,probs = 0.975),
  s.010<-  apply(x,y,quantile,probs = 0.10),
  s.090<-  apply(x,y,quantile,probs = 0.90))}

t.1=sumr(occ1$psi)
t.2=sumr(occ2$psi)
a=list()

# 
occ.ls <- list(list(t.1,"cond1"),list(t.2,"baseline"))

for (i in 1:nspecies){
  x1 <- list()
  for (c in 1:length(occ.ls)){
    x.1  <- occ.ls[[c]][[1]]
    x.2  <- occ.ls[[c]][[2]]
    x1[[c]]=data.frame(
      Mean=     x.1[[1]][i,],
      LB_95  =  x.1[[2]][i,], 
      UB_95 =   x.1[[3]][i,],
      LB_80 =   x.1[[4]][i,], 
      UB_80 =   x.1[[5]][i,],
      Year = seq(1,10,1),
      Cond = rep(x.2,10)
    )
  }
  x1  <- data.table::rbindlist(x1)
  a[[i]]=plot_div(x1,title=NULL,subtitle = "Trend", color="#3399FF",
                  x.axis.texts = element_text(size=8),x.title.texts =element_text(size=8),
                  y.axis.texts = element_text(size=8),group=T)
}
a

# maps
# percent or 0 beta


beta.s<- lapply(b.pair.x,beta_check)
beta.sm<- lapply(b.mult.x,beta_check,M=T)
list(beta.s,beta.sm)


l.x=beta.sm
bta<- list(list(),list(),list())

for(l in 1:length(l.x)){
  for (i in 1:3){
    l.x[[l]][[i]]$Cond <- paste0(l)
    l.x[[l]][[i]]$Year <- 1:10
  }}

ax <- data.table::rbindlist(l.x[[1]])

plot_div(l.x[[1]][[1]],title=NULL,subtitle = "Trend", color="#3399FF",
         x.axis.texts = element_text(size=8),x.title.texts =element_text(size=8),
         y.axis.texts = element_text(size=8),group=T)
###############################################################################
# pairs 

t1 = psi[,,,1]
t2 = psi[,,,2]
b=dim(psi)[c(2,3,1)]
comb.arr <- array(dim=c(2,b[1],b[2],b[3]))
for ( s in 1:b[1]){
  for ( i in 1:b[2]){
    for ( nits in 1:b[3]){
      comb.arr[1,s,i,nits] <- t1[nits,s,i]
      comb.arr[2,s,i,nits] <- t2[nits,s,i]
    }}}
comb.arr[,,1,1]

####################################################################
beta_div <- function(x,occ=F){
  
  if(occ==T){
    b.mult.x <- apply(aperm( x,c(3,2,1,4)),c(3,4),
                      function(y){betapart::beta.multi.abund(y)})
    
    b.pair.x <- apply(aperm( x,c(4,2,1,3)),c(3,4),
                      function(y){betapart::beta.pair.abund(y)})
  } else {
    
    b.mult.x <- apply(aperm( x,c(3,2,1,4)),c(3,4),
                      function(y){betapart::beta.multi(y)})
    
    b.pair.x <- apply(aperm( x,c(4,2,1,3)),c(3,4),
                      function(y){betapart::beta.pair(y)})
  }
list(b.mult.x,b.pair.x)

}


#
################################################################################
beta_check<- function(t,M=F){
  
  betas<-apply(t,c(2),function(x){x})
  
  if(M!=T){
    x1=lapply(betas, function(x) unlist( lapply(x,function(x) as.matrix(x[[1]])[f.year,s.year])))
    x2=lapply(betas, function(x) unlist( lapply(x,function(x) as.matrix(x[[2]])[f.year,s.year])))
    x3=lapply(betas, function(x) unlist( lapply(x,function(x) as.matrix(x[[3]])[f.year,s.year])))
  }else {
    
    x1=lapply(betas, function(x) unlist( lapply(x,function(x) x[[1]])))
    x2=lapply(betas, function(x) unlist( lapply(x,function(x) x[[2]])))
    x3=lapply(betas, function(x) unlist( lapply(x,function(x) x[[3]]))) 
    
  }
  
  beta.s <- lapply(list(x1,x2,x3),function(x){x2  <- lapply(x,function(x)
    data.frame(Mean=mean(x),  
               LB_95=quantile(x,probs = 0.025),
               UB_95=quantile(x,probs = 0.975),
               LB_80=quantile(x,probs = 0.10),
               UB_80=quantile(x,probs = 0.90)))
  
  
  x3 <- data.table::rbindlist(x2)
  
  if (M!=T){
    x3$gr <- names(x2)
    coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
    left_join(x3,coord)}
  else{
    
    x3}})
  
  beta.s
}
###########################################################################
s.year=1
f.year <- dim(occ.s$psi)[[4]]
reps <- length(x)

occ.s <- readRDS("occ.1.rds")
t=b[[2]]

beta.s<- beta_check(b[[2]],M=F)
beta.sm<- lapply(b.mult.x,beta_check,M=T)

##########################################################################
############### Parameter effects across England #########################
##########################################################################
x=b[[2]][[1]]
x1=beta.s[[1]]
ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = beta.s[[3]], 
                               aes(x = E, y = N, fill =UB_95))+
  scale_fill_continuous(type = "viridis", name = "Similarity")


b <- readRDS("beta.rds")
