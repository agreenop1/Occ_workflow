################################################################################
################################################################################
##################### Autocorrelation ##########################################
library(ggplot2)
d1 <- readRDS("jas_out/summary_c.rds")


# summarize output
summ <- function(x){
  m1 <- apply(x,c(2),mean,na.rm=T)
  m2 <- apply(x,c(2),quantile,probs = 0.025,na.rm=T)
  m3 <- apply(x,c(2),quantile,probs = 0.975,na.rm=T)
 
  trend <- data.frame(mean=m1,lb=m2,ub=m3)
  trend$year <- seq(10,680,10)
  trend
}

do <- summ(d1$spatial$observed)
do$Trend <- "observed"
de <- summ(d1$spatial$expected)
de$Trend <- "expected"
  
dt <- rbind(do,de)

c1 <- ggplot(data=dt)+geom_point(aes(x=year,y=mean,color=Trend))+geom_line(aes(x=year,y=mean,color=Trend))+
  geom_ribbon(aes(x=year,ymin=lb,ymax=ub,fill=Trend),alpha=0.2)+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Correlation")+xlab("Distance (km)")

ggsave('check_plot/auto_cor_plots.png',c1,width=8,height=6)


t1 <- d1[[2]][[1]]
f <- apply(t1,2,function(x){ t2<- do.call(rbind,x)
                       colMeans( matrix(ncol=11,nrow=2278, unlist(t2[,1]),byrow = T))})

# summarize output
summt <- function(x){
  m1 <- apply(x,c(1),mean,na.rm=T)
  m2 <- apply(x,c(1),quantile,probs = 0.025,na.rm=T)
  m3 <- apply(x,c(1),quantile,probs = 0.975,na.rm=T)
  
  trend <- data.frame(mean=m1,lb=m2,ub=m3)

  trend
}
f <- summt(f)
c1 <- ggplot(data=f)+geom_point(aes(x=1:11,y=mean,color="red"))+
  geom_errorbar(aes(x=1:11,ymin=lb,ymax=ub,color="red"))+geom_hline(yintercept = 0,linetype="dashed")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  )+ylab("Correlation")+xlab("Lag")
ggsave('check_plot/temp_cor_plots.png',c1,width=8,height=6)
