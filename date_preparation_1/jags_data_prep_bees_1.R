#########################################################################################
# 1. Script and functions to match occupancy in the BRC format to environmental         #
# covariates to be used in dynamic species occupancy models                             #
#########################################################################################
# Load packages
library(sparta)
library(dplyr)
library(BRCmap)
library(stringr)
source("brcmap_f.R")
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
years_n <- seq(1994,2016,2)
UK <-  readRDS("UK_map.rds")
plots <- F
# closure period 1 = yearly and 2 = biennial 
# Lag determines whether the covariates for t (lag = T) or t-1 are used to predict persistence

closure.period=2

# If lag = T then occupancy (z) for time period t is dependent on z in t-1 & persistence driven by covariates in t-1.
# Here t is either equal to a year or every 2 years dependent on closure period.
# We need occurence records spanning 1st time period we have covariates (t1) to the last time period + t (ti+t),
# as the covariates are being used to predict z in the following time period t+1.

# If lag = F then z for time period t is dependent on z in t-1 & persistence driven by covariates in t.
# We need occupancies spanning 1st time period we have covariates for - t (t1-t) 
# to the last time we have covariates for (ti), as we use z in t-1 and covariates in t 
# to predict z in t.

lag = T

# group
occ.dat = readRDS("occurrence_recs/bee.rds")
group.name="bees"

# sites at the 5km resolution  
occ.dat$grid <- nchar(as.character(occ.dat$TO_GRIDREF))

occ.dat <- occ.dat[occ.dat$grid>5,]

# years
start_year =1994; end_year = 2016

# reformat to 5km grid
occ.dat$site_5km <- reformat_gr(occ.dat$TO_GRIDREF, 
                                prec_out = 5000) # in your case you want prec_out = 5000
# threshold for species inclusion 
# all = all species; com.20 = common species (.n selects number e.g. com.20 = 20 most common species); 
# rar.20 = rare species (require min of 50 obs) ; type= list(list of species names-for specific species)

type="all"

n.sample=F # sample observations - prop

obs.n=499 # min threshold for species inclusion where type = all

################################################################
# covars
chem_cov <- read.csv("raw_FERA/species_RQ.csv") # chem data

# env covariates
env_cov <- readRDS("env.cov.list.rds")

# sum all risk quotients 
rq_sum <-  chem_cov %>% group_by(year,gr,E,N)%>% 
           summarise(predators_RQ=sum(predators_RQ,na.rm=T),
           predators_pollinators_RQ=sum(predators_pollinators_RQ,na.rm=T),
           pollinators_RQ=sum(pollinators_RQ,na.rm=T))

# get the mean of the risk quotient
mean_RQ <- rq_sum %>% group_by(gr) %>% summarise(rq_M_pol= mean(pollinators_RQ),
                                                 rq_M_pre= mean(predators_RQ),
                                                 rq_M_ppl= mean(predators_pollinators_RQ))

# separate risk quotient into spacial and temporal variable
RQ <- left_join(rq_sum,mean_RQ)

RQ$rq_A_pol <- RQ$pollinators_RQ - RQ$rq_M_pol
RQ$rq_A_pre <- RQ$predators_RQ - RQ$rq_M_pre
RQ$rq_A_ppl <- RQ$predators_pollinators_RQ - RQ$rq_M_ppl

# zero application - used in simulations
zero_app <- data.frame(RQ[1:4],rqpol=0-RQ$rq_M_pol,
                       rqpre=0-RQ$rq_M_pre,
                       rqppl=0-RQ$rq_M_ppl)

# covs to be analysed - these are passed to var_prep function (see below):
# site is currently set as "gr" can be changed in the function
# need to be in long format with a year and gr column then the variable of interest as var.x below

cov_assess <- c(list(temp_anom=env_cov$temp_anom,
                     mean_temp=env_cov$mean_temp,
                     semi=env_cov$semi, 
                     agri=env_cov$agri),
                     list(RQsum_A=RQ),
                     list(RQsum_M=RQ)) # named list of covs

# do variables need to be centered or not
scale_bin <- c(F,T,T,T,F,T)

var.x <- c("mean.x","mean.x","mean.x","mean.x",
           "rq_A_pol","rq_M_pol") # variable of interest CHECK!


#################################################################  
# derive site list  based on common sites in the covariate data and occupancy data
occ.site <- occ.dat$site_5km
sites <-function(){list(unique(occ.site), unique(chem_cov$gr),unique(env_cov[[1]]$gr)) }
sites_5km <- Reduce(intersect, sites()) # common list of sites

# FUNCTION #######################################################
# covariates in wide format
# checks if the order of sites in the covariate data match there numerical order in occ data
# if not then reorder covariate data to ensure it matches
# keep.id = keep column with site id (needed for checks)
# center = center variables?

cov_year <- seq(start_year,end_year,2) # years

var_prep <- function(x,var.x,var.names=NULL,site="gr",keep.id=F,center=F){
  
  cat(var.names,"\n")
  x <- x[x$gr%in%site_id$site_name,] # make sure only sites in occ are used
  x <- x[x$year%in%cov_year,]
  x <- x[c("year",site,var.x)]
  
  # get every site x year combination to check for missing/0 values
  yxs <- expand.grid(sites_5km,cov_year) # all combos
  colnames(yxs)[1:2] <- c("gr","year")
  yxs$comb <- paste0(yxs$gr,"_",yxs$year)
  
  miss_comb   <- yxs[!yxs$comb%in%paste0(x$gr,"_",x$year),] # list of missing y x s
  
  # if y x s combos is greater than 0 
  if(nrow(miss_comb)>0){
    fill <- readline(prompt = "Missing values, fill with 0? (yes/no)")
    
    if(fill=="yes"){
      miss_comb[,var.x] <- NA # create new cols
      miss_comb$comb <- NULL
      x <- rbind(x,miss_comb[-5]) # bind the missing sites with NA
      }else{
      cat("Data contains missing (NA) values!")
      }
  }
  
  if(center){x[var.x] <- scale(x[var.x],scale = F)} # center var if T
  
  x1 <- pivot_wider(data=x,id_cols=site,names_from=year,values_from = var.x) # make wide
  x1 <- x1[c("gr",cov_year)] # keep years we need
  
  # should na be filled with 0?!
    if (any(is.na(x1[paste0(cov_year)]))){
      if(fill=="yes"){cat("  Filling NA with 0!","\n")
      x1[is.na(x1)] <- 0 }else {cat("Not filling NA","\n")}}
  
  # output from fill
  x= x1
  
  # make sure sites are in the correct order as occ
  if(all(as.character(x$gr)==site_id$site_name)){
    cat("  Site order is correct","\n")
    x
  }else{
    cat("  Site order is INCORRECT. Rearranging covariate data to ensure it matches site numbers.","\n")
    x <- x[match(site_id$site_name,x$gr),] # if order is incorrect put them in the correct order
    
  }
  if(!keep.id){x[,-1]}else{x}
}


##########################################################################
# main code for outputs #
##########################################################################
# is a sample selected 
if(is.numeric(n.sample)){ occ.dat <- occ.dat[sample(1:nrow(occ.dat),nrow(occ.dat)*n.sample),]}


#  if function to ensure the correct occurrence records are selected across years
years <- if(lag&closure.period==1){start_year:(closure.period+end_year)
         }else if(lag&closure.period==2) {start_year:(closure.period+end_year+1)
         }else if(!lag&closure.period==1) {(start_year-closure.period):(end_year)
         }else {(start_year-closure.period):(end_year+1)}


##########################################################################
# subset occ.dat on common sites
occ.dat <- occ.dat[occ.dat$site_5km%in%sites_5km,]

#remove those that don't have day precision
occ.dat$TO_STARTDATE <- as.Date(occ.dat$TO_STARTDATE, format = "%Y-%m-%d")

occ.dat$TO_ENDDATE <- as.Date(occ.dat$TO_ENDDATE, format = "%Y-%m-%d")

occ.dat$DT_ID <- occ.dat$TO_ENDDATE - occ.dat$TO_STARTDATE

occ.dat <- occ.dat[occ.dat$DT_ID<1,]


# remove years were not interested in
occ.dat <- occ.dat[occ.dat$YEAR%in%years, ] 

###########################################################################
# closure period 
# if closure period = 2 observations are grouped together at 2-year intervals
if(closure.period==2){
  
  n_periods <- rep(1:(length(years)/closure.period),each=2) # Create vector half the length of years
  clsr_per <- data.frame(n_periods,years) #  assign new temporal ID to each year
  colnames(clsr_per)[1:2] <- c("closure_per","YEAR")
  occ.dat <- left_join(occ.dat,clsr_per) #  Match in the main occupancy data year to new closure period
  
  print (distinct(occ.dat[c("YEAR","closure_per")]))}else{

  # yearly closure period 
occ.dat$closure_per <- occ.dat$YEAR-(min(occ.dat$YEAR-1))

}

################ NEED TO RUN TO GET CORRECT LIST LENGTH #################
  #  set up data for JAGS format
  # format data
  occ <- formatOccData(taxa = occ.dat$CONCEPT,
                       site = occ.dat$TO_GRIDREF,
                       survey =  occ.dat$TO_STARTDATE,
                       closure_period = occ.dat$closure_per)
  
# reformat to 5km grid
occ[[2]]$site_5km <- reformat_gr(occ[[2]]$site, 
                                 prec_out = 5000) # 5km2

#  create separate date frames for the occupancy info & visit info
occup <- occ[[1]]
visit <- occ[[2]]
occup[occup==T] <- 1

#########################################################################
# remove sites only visited in one TP
site_vis <- distinct(visit[c("site","TP")])
site.v.n1  <- data.frame (table(site_vis$site))
site.v.n2  <- site.v.n1[site.v.n1$Freq>1,]

# 
visit <- visit[visit$site%in%site.v.n2$Var1,]
occup <-  occup[occup$visit%in%visit$visit,]
  
# species selection #####################################################
#  filter out species with fewer than obs.n observations


  occ.n <- data.frame(sp=colnames(occup[-1]),tot= colSums(occup[-1]))
  occ.n <- occ.n[occ.n$tot>obs.n,]
  occup <- occup[c("visit",paste(occ.n$sp))]
  f.name=paste0("all.",obs.n)

colSums(occup[-1])

#########################################################################
# this is run again as sites are removed
# derive site list  based on common sites in the covariate data and occupancy data
occ.site <- unique(visit$site_5k)
sites_5km <- Reduce(intersect, sites())

#########################################################################
#  set site levels numerically &  ensure covariate data is in the same order (based on sites)

visit$site_5km.n <-  as.numeric(droplevels(as.factor(visit$site_5km))) # assign each site a numerical value

# get each unique site value, actual site name and order based on the numerical value of site in occ data

site_id <- distinct(visit[5:6])
colnames(site_id) <- c("site_name","site_num")
site_id <- site_id[order(site_id$site_num), ] 

#########################################################################
# observation model list length effects 
visit$SHORT <- 0
visit$SHORT[between(visit$L,2,3)] <- 1

visit$LONG <- 0
visit$LONG[visit$L>3] <- 1


visit$Year <- as.numeric(str_sub(visit$visit,start=7,end = 10 ))-1993
nYear <- length(unique(visit$Year))
#########################################################################
# make sure occupancy and visit data match
if(all(occup$visit==visit$visit)){"Occupancy and visit data match"}

# zobs
obs <- cbind(occup,visit)
 
obs.l <- pivot_longer(obs,cols =colnames(occup[-1]),names_to = "CONCEPT" )
zobs <- acast(obs.l,CONCEPT~site_5km~TP,value.var = "value",fun.aggregate=sum, fill = -9999)
zobs[zobs>1] <- 1
zobs[zobs==-9999] <- NA
zobs[zobs==0] <- NA # stops conflicts in nimble

########################################################################
# covariates with explicit site ID need for checks
var.names = names(cov_assess)
covars_id <- list()

for (i in 1:length(var.names)){
covars_id[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=T,center=F)
names(covars_id)[i] <- var.names[[i]]
}


# covariates used in simulations
z1 <- var_prep(zero_app,var.x="rqpol", site="gr", keep.id=T,center=F) # zero application wide format
z0 <- var_prep(cov_assess$RQsum_A,var.x="rq_A_pol", site="gr", keep.id=T,center=F) # actual application rq temporal
saveRDS(list(zero_application=z1,actual=z0),"zero_app.rds") 
write.csv(covars_id$temp_anom,"site_id.csv")
saveRDS(covars_id,"covars_wide.rds")

# double check all covs are in the right order
for (i in 1:length(covars_id)){
  if(!all(covars_id[[i]]$gr==site_id$site_name)){stop("covs are not in the right order")}
}

# covariates without site id
covars <- list()

for (i in 1:length(var.names)){
  covars[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=F,center=scale_bin[[i]])
  names(covars)[i] <- var.names[[i]]
}
#########################################################################
if(plots){
  # species plots
  coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])
  hover <- read.csv("hoverfly_names.csv")
  
  keep_names <- hover$species[hover$tribe%in%unique(hover$tribe)] # Bacchini Eristalini
  
  occ_record <- cbind(occ_tot=rowSums( occup[keep_names]),visit)
  occ_map <- occ_record %>% group_by(TP,site_5km) %>% summarise(occ_sum=occ_tot)
  colnames(occ_map)[2] <- 'gr'
  
  occ_map <- left_join(occ_map,coord)
  
  
  out <- list()
  
  for(i in 1:12){
    
    sy <- occ_map[occ_map$TP==i,]
    
    sy <- sy[sy$occ_sum>0,]
    
    
    out[[i]] <- ggplot() +
      geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
      ylim(0, 700000)  + geom_tile(data = sy, 
                                   aes(x = E, y = N, fill =occ_sum))+
      scale_fill_continuous(type = "viridis", name = years_n[[i]])+ theme(panel.grid.major = element_blank(), 
                                                                          panel.grid.minor = element_blank(),
                                                                          panel.background = element_blank())
  }
  
  out
}
########################################################################
# covariates with explicit site ID need for checks
var.names = names(cov_assess)
covars_id <- list()

for (i in 1:length(var.names)){
  covars_id[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=T,center=F)
  names(covars_id)[i] <- var.names[[i]]
}


# covariates used in simulations
z1 <- var_prep(zero_app,var.x="rqpol", site="gr", keep.id=T,center=F) # zero application wide format
z0 <- var_prep(cov_assess$RQsum_A,var.x="rq_A_pol", site="gr", keep.id=T,center=F) # actual application rq temporal
saveRDS(list(zero_application=z1,actual=z0),"zero_app_hovpol.rds") 
write.csv(covars_id$temp_anom,"site_id.csv")
saveRDS(covars_id,"covars_wide_hovpol.rds")

# double check all covs are in the right order
for (i in 1:length(covars_id)){
  if(!all(covars_id[[i]]$gr==site_id$site_name)){stop("covs are not in the right order")}
}

# covariates without site id
covars <- list()

for (i in 1:length(var.names)){
  covars[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=F,center=scale_bin[[i]])
  names(covars)[i] <- var.names[[i]]
}

#########################################################################
#################### Risk of Bias Assessment ############################
#########################################################################
# recording efforts
sxc <- distinct( visit[c("TP","site_5km")])
n.svisits <-  as.data.frame( table(sxc$site_5km)) # number of time period sites visited in
colnames(n.svisits)[1] <- "site_5km"

# region observation 
region <-distinct( read.csv("osr_pesticide_5km_1994_2010_v3.csv")[c("ref_5km","region")]) 
region$region <- tolower(region$region)
region <-distinct(region) 
colnames(region)[1] <- "site_5km"

n.svisits.r <- left_join(n.svisits,region)
visits_reg <- left_join(visit,region)

# frequency of observations at the uk and regional scale
uk_freq <- do.call(rbind, lapply(2:13,function(x) data.frame(percent=(sum(n.svisits.r$Freq==x )/sum(n.svisits.r$Freq>0 ))*100,TP_visited=x)))

region_freq <- do.call(rbind, lapply(2:13,function(x){
                                      region_freq <- n.svisits.r %>% group_by(region) %>% summarise(percent=(sum(Freq==x)/sum(Freq>0)*100))
                                      region_freq$TP_visited <- x
                                      region_freq$region <- as.factor(region_freq$region )
                                      region_freq
                                      }))
                        
# plots of sampling frequency at different scales
# uk
uk_plot_f <- ggplot()+ geom_line(data=uk_freq,aes(y=percent,x=TP_visited))+
          geom_point(data=uk_freq,aes(y=percent,x=TP_visited))

# regional
region_plot_f <- ggplot()+ geom_line(data=region_freq,aes(y=percent,x=TP_visited,colour=region))+
        geom_point(data=region_freq,aes(y=percent,x=TP_visited,colour=region))

############################################
# number of observations at different scales
region_v <- as.data.frame(table(visits_reg$region,visits_reg$TP))
site_v <- as.data.frame(table(visits_reg$site_5km))
site_f <- distinct(visit[c("site_5km","TP")])
site_f <- as.data.frame(table(site_f$site_5km))

# set column names
colnames(region_v) <- c("region","Closure","Number_of_visits")
colnames(site_v) <- c("gr","number_of_visits")
colnames(site_f) <- c("gr","number_of_closure")

# create different regional parameters
region_v_c <- region_v %>% group_by(Closure) %>% summarise(total=sum(Number_of_visits))
region_v <- left_join(region_v ,region_v_c )
region_v$Closure <- as.numeric(region_v$Closure)

# regional plots of the total number of visits
region_plot_v <- ggplot()+ geom_line(data=region_v,aes(y=Number_of_visits ,x=Closure,colour=region))+
  geom_point(data=region_v,aes(y=Number_of_visits ,x=Closure,colour=region))

region_plot_v_stacked <-   ggplot( ) + 
  geom_bar(data=region_v,aes( y=Number_of_visits, x=Closure,fill=region),position="stack", stat="identity")

##############################
# map of total visits per site
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4]) # read in coordinates for plots 
site_v  <- left_join(site_v,coord) 
site_f  <- left_join(site_f,coord) 

# total number of visits by 5km cell
ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = site_v , 
                               aes(x = E, y = N, fill =number_of_visits))+
  scale_fill_continuous(type = "viridis", name = 'Total visits')+ theme(panel.grid.major = element_blank(), 
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.background = element_blank())
# total number of closure periods Visited in by 5km cell
ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = site_f , 
                               aes(x = E, y = N, fill =number_of_closure))+
  scale_fill_continuous(type = "viridis", name = 'Total number of closure periods \n visited in')+ 
  theme(panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank())

#################################
# look at sampling across species 
occ.dat1 <- occ.dat[occ.dat$CONCEPT %in% colnames(occup[-1]),]
occ.dat1 <- left_join(occ.dat1,region)
sp_re_closure <- distinct( occ.dat1[c("CONCEPT","region","closure_per")])
sp_si_closure <- distinct( occ.dat1[c("CONCEPT","site_5km","closure_per")])
sp_closure <- distinct( occ.dat1[c("CONCEPT","closure_per")])

# table different variables
n_sp_cl <- as.data.frame( table(sp_closure$closure_per))

n_sp_ev <- as.data.frame( table(sp_si_closure$CONCEPT,sp_si_closure$closure_per))
n_sp_ev_w<- pivot_wider(n_sp_ev, names_from=Var2, values_from=Freq )            

n_spre_cl <- as.data.frame( table(sp_re_closure$region,sp_re_closure$closure_per))

# name columns
colnames(n_sp_cl) <- c("Closure","n_sp_observed")
n_sp_cl$Closure <- as.numeric(n_sp_cl$Closure)
colnames(n_spre_cl) <- c("region","Closure","n_sp_observed")
n_spre_cl$Closure <- as.numeric(n_spre_cl$Closure)

ggplot()+ geom_line(data=n_sp_cl,aes(y=n_sp_observed,x=Closure))+
  geom_point(data=n_sp_cl,aes(y=n_sp_observed,x=Closure))

ggplot()+ geom_line(data=n_spre_cl,aes(y=n_sp_observed,x=Closure,color=region))+
  geom_point(data=n_spre_cl,aes(y=n_sp_observed,x=Closure,color=region))


#########################################################################
########################### covariate plots #############################  
#########################################################################
if(plots){
  # plots 
  
  
  id_plots <- lapply(covars_id[1:6],function(x){left_join(x,coord)}) # join covariates with coordinates
  
  # time x space plots for all covariates
  env_plots <- lapply (1:6,function(x1){
    
    x1 <- id_plots[[x1]]
    out <- list()
    
    for(i in 1:12){
      
      sy <- x1[c("E","N",years_n[[i]])]
      
      sy <- sy[sample(1:2000,500),]
      colnames(sy)[1:3] <- c("E","N","Year")
      
      out[[i]] <- ggplot() +
        geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
        ylim(0, 700000)  + geom_tile(data = sy, 
                                     aes(x = E, y = N, fill =Year))+
        scale_fill_continuous(type = "viridis", name = years_n[[i]])+ theme(panel.grid.major = element_blank(), 
                                                                            panel.grid.minor = element_blank(),
                                                                            panel.background = element_blank())
    }
    
    out
  })
  
  names(env_plots) <- names(id_plots)
  
  p1 <- env_plots
  ggsave("pmap.png",plot=p1, dpi=600, dev='tiff',   height=4.5, width=4.5, units="in") 
  
  # FUNCTION ##############################################################
  #  Check for temporal correlation in parameters at sites through years
  pair <- combn( names(id_plots),2)
  
  n_sites  <- nrow(site_id)  
  cor_out <- data.frame(matrix(nrow=n_sites,ncol=ncol(pair)))
  
  for(j in 1:ncol(pair)){
    
    cname <- pair[,j]
    colnames(cor_out)[j] <- paste(cname[1],"vs",cname[2])
    
    x1 = id_plots[[cname[1]]];y1 = id_plots[[cname[2]]]
    x1 = x1[paste0(years_n)]; y1 = y1[paste0(years_n)]
    
    
    
    
    for(i in 1:n_sites){
      cor_out[i,j] <-  cor(t(x1[i,]),t(y1[i,]))
    }  
  }
  
  # correlation plots - temporal
  rt <- ggplot(cor_out,aes(x=`temp_anom vs RQsum_A`)) + geom_histogram()+
    xlab("Temperature anomaly vs RQ anomaly (r)")
  
  ar <- ggplot(cor_out,aes(x=`agri vs RQsum_A`)) + geom_histogram()+
    xlab("Agriculture vs RQ anomaly (r)")
  
  ta <- ggplot(cor_out,aes(x=`temp_anom vs agri`)) + geom_histogram()+
    xlab("Agriculture vs Temperature anomaly (r)")
  
  st <- ggplot(cor_out,aes(x=`temp_anom vs semi`)) + geom_histogram()+
    xlab("Semi-natural vs Temperature anomaly (r)")
  
  sr <- ggplot(cor_out,aes(x=`semi vs RQsum_A`)) + geom_histogram()+
    xlab("Semi-natural vs RQ anomaly (r)")
  
  sa <- ggplot(cor_out,aes(x=`semi vs agri`)) + geom_histogram()+
    xlab("Semi-natural vs Agriculture (r)")
  
  temr_plots <- ggarrange(plotlist=list(rt,ar,ta,st,sr,sa),nrow=3,ncol=2)
  
  ggsave("check_plot/temporal_plot.png",temr_plots, width = 7,height=8)
  
  # FUNCTION ##############################################################
  # check for correlation between parameters
  # get combos of all variables
  pair <- combn(names(covars_id[1:6]),2)
  x = covars_id
  # calculate correlations      
  cors_mat  <- apply(pair,2,function(cname){ 
    x1 = x[[cname[1]]];y1 = x[[cname[2]]]
    
    x2  <- pivot_longer(x1,colnames(x1[-1]))
    y2  <- pivot_longer(y1,colnames(y1[-1]))
    colnames(y2)[3] <- "value1"
    cors <-  matrix(ncol =4,nrow=ncol(x1[-1]))
    c3 <- left_join(x2,y2)
    pcor <-  ggplot(data=c3)+geom_point(aes(x=value,y=value1,color=as.numeric(name)))+ 
      xlab(cname[[1]]) + ylab(cname[[2]])
    pcor$labels$colour <- "Year"
    pcor
    
    for( i in 1:ncol(x1[-1])){
      
      cors[i,1] <- i
      cors[i,2]<-  round(cor(x1[i+1],y1[i+1]),2)
      cors[i,3] <-cname[1]
      cors[i,4] <-cname[2]
    }
    
    
    list(cors,pcor)
  })
  
  
  #  correlation plot + table of correlations
  lapply(cors_mat,function(x){
    cmx <-as.data.frame(x[[1]])
    colnames(cmx) <- c("Year","Cor","Var 1", "Var 2")
    cmx$Year <- cov_year
    plt <- x[[2]]
    plts <- ggarrange(plt,ggtexttable(as.data.frame(cmx)),nrow=2)
    ggsave(paste0("check_plot/",cmx[1,3],"_",cmx[1,4],".png")
           ,plot=plts,width = 8,height = 8)
  })
  
  # mean temporal trends
  
  time_trends <- lapply(covars_id,function(x1){
    cv.trend  <-  pivot_longer(x1,colnames(x1[-1])) %>% 
      group_by(name) %>% 
      summarise(M=mean(value),SD=sd(value))
    
    ggplot(data=cv.trend)+
      geom_line(aes(x=name,y=M,group=1))+
      geom_point(aes(x=name,y=M,group=1))+
      geom_errorbar(aes(x=name,ymin=M-SD,ymax=M+SD))})
  
  ggsave("check_plot/RQts.png",plot=time_trends$RQsum,width = 8,height = 2)
  ggsave("check_plot/AGts.png",plot=time_trends$agri,width = 8,height = 2)
  ggsave("check_plot/SEts.png",plot=time_trends$semi,width = 8,height = 2)
  
  ##########################################################################
  # variable trends over time
  time_trends <- lapply(covars_id,function(x){
    
    # remove tile column
    cx <- x[-1]
    cx_out <- data.frame(matrix(ncol = 3,nrow = 12))
    
    # get mean and 95 credible intervals
    cx_out[,1] <- apply(cx,2,mean)
    cx_out[,2] <- apply(cx,2,quantile,probs = 0.025)
    cx_out[,3] <- apply(cx,2,quantile,probs = 0.975)
    cx_out
  })
  
  
  
  tpt <- ggplot(time_trends[[1]],aes(x=seq(1994,2016,2),y=X1))+geom_line()+geom_errorbar(aes(ymin=X2, ymax=X3)) +xlab("Years")+
    ylab("Temperature anomaly")
  ggsave("check_plot/time_t.png",plot=tpt,width = 8,height = 2)
  rqt <- ggplot(time_trends[[5]],aes(x=seq(1994,2016,2),y=X1))+geom_line()+geom_errorbar(aes(ymin=X2, ymax=X3)) +xlab("Years")+
    ylab("Risk quotient anomaly")
  ggsave("check_plot/rq_t.png",plot=rqt,width = 8,height = 2)
  
  
  
}
#########################################################################
# order check
if(!all(rownames(zobs[1,,])==site_id$site_name)) { cat("check sites are in the correct order")}

for (i in 1:length(covars_id)){
  if(!all(covars_id[[i]]$gr==rownames(zobs[1,,])))  {cat("check sites match covariate sites","\n")}
}

if(!all(dimnames(zobs)[[1]]==colnames(occup[-1]))){cat("check species names are in the correct order")}

#########################################################################
# region observation 
region <-distinct( read.csv("osr_pesticide_5km_1994_2010_v3.csv")[c("ref_5km","region")]) 
region$region <- tolower(region$region)
region <-distinct(region) 
region$region.n <- as.numeric(as.factor(region$region))
colnames(region)[1] <- "site_5km"
visit <- left_join(visit,region[c("site_5km","region")])

# region ecological
region.cov <- left_join(data.frame(site_5km=rownames(zobs[1,,])),region)
nregion <- length(unique(region.cov$region.n))
# order check
if(!all(rownames(zobs[1,,])==region.cov$site_5km)) { cat("check sites are in the correct order")}

#########################################################################
# save outputs
cat("Minimum observations",obs.n,"\n")
cat(nrow(site_id),"sites","\n")
cat(ncol(occup[-1]),"species","\n")

cat(file.name <- paste0("Model_data/data_",group.name,"_",f.name,"_",start_year,".",end_year,".rds"),"\n")



saveRDS(list(occup,visit,zobs,covars,closure.period,nYear=nYear,region=region.cov$region.n,nregion=nregion),file.name)
jas_data <- list(occup,visit,zobs,covars,closure.period,nYear=nYear,region=region.cov$region.n,nregion=nregion)
jas_data$region
jas_data$nregion
