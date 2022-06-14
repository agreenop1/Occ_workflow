#########################################################################################
# 1. Script and functions to match occupancy in the BRC format to environmental         #
# covariates to be used in dynamic species occupancy models                             #
#########################################################################################
# Load packages
library(sparta)
library(dplyr)
library(BRCmap)
library(dplyr)
library(reshape2)
library(tidyr)


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
start_year =1994; end_year = 2010

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
chems <- split(chem_cov,f=chem_cov$MOA)
env_cov <- readRDS("env.cov.list.rds")

# sum all risk quotients 
rq_sum <-  chem_cov %>% group_by(year,gr,E,N)%>% summarise(predators_RQ=sum(predators_RQ,na.rm=T),
                                                           predators_pollinators_RQ=sum(predators_pollinators_RQ,na.rm=T),
                                                           pollinators_RQ=sum(pollinators_RQ,na.rm=T))


# covs to be analysed - these are passed to var_prep function (see below):
# site is currently set as "gr" can be changed in the function
# need to be in long format with a year and gr column then the variable of interest as var.x below

cov_assess <- c(list(temp_anom=env_cov$temp_anom,
                     mean_temp=env_cov$mean_temp,
                     semi=env_cov$semi, 
                     agri=env_cov$agri),
                     list(RQsum=rq_sum)) # named list of covs

var.x <- c("mean.x","mean.x","mean.x","mean.x","pollinators_RQ") # variable of interest CHECK!

#################################################################  
# derive site list  based on common sites in the covariate data and occupancy data
sites <-function(){list(unique(occ.dat$site_5km), unique(chem_cov$gr),unique(env_cov[[1]]$gr)) }
sites_5km <- Reduce(intersect, sites())

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


################ NEED TO RUN TO GET CORRECT LIST LENGTH #################
#  set up data for JAGS format

# if closure period = 2 observations are grouped together at 2-year intervals
if(closure.period==2){
  
  n_periods <- rep(1:(length(years)/closure.period),each=2) # Create vector half the length of years
  clsr_per <- data.frame(n_periods,years) #  assign new temporal ID to each year
  colnames(clsr_per)[1:2] <- c("closure_per","YEAR")
  occ.dat <- left_join(occ.dat,clsr_per) #  Match in the main occupancy data year to new closure period
  
  print (distinct(occ.dat[c("YEAR","closure_per")]))
  

  # remove sites only visited in one TP
  site_vis <- distinct(occ.dat[c("TO_GRIDREF","closure_per")])
  site.v.n1  <- data.frame (table(site_vis$TO_GRIDREF))
  site.v.n2  <- site.v.n1[site.v.n1$Freq>1,]
 
  # 
  occ.dat <- occ.dat[occ.dat$TO_GRIDREF%in%site.v.n2$Var1,]
  
  # format data
  occ <- formatOccData(taxa = occ.dat$CONCEPT,
                       site = occ.dat$TO_GRIDREF,
                       survey =  occ.dat$TO_STARTDATE,
                       closure_period = occ.dat$closure_per)
}else{
  
  # remove sites only visited in one TP
  site_vis <- distinct(occ.dat[c("TO_GRIDREF","YEAR")])
  site.v.n1  <- data.frame (table(site_vis$TO_GRIDREF))
  site.v.n2  <- site.v.n1[site.v.n1$Freq>1,]
  
  # 
  occ.dat <- occ.dat[occ.dat$TO_GRIDREF%in%site.v.n2$Var1,]
  
  
  #   use year as the closure period
  # format data
  occ <- formatOccData(taxa = occ.dat$CONCEPT,
                       site = occ.dat$TO_GRIDREF,
                       survey =  occ.dat$TO_STARTDATE)
}

# reformat to 5km grid
occ[[2]]$site_5km <- reformat_gr(occ[[2]]$site, 
                                 prec_out = 5000) # 5km2

#  create separate date frames for the occupancy info & visit info
occup.f <- occ[[1]]
visit.f <- occ[[2]]

# species selection ##########################################################
#  filter out species with fewer than obs.n observations
if(type=="all"){
  cat("All species selected","\n")
  obs.n=obs.n
  occ.dat <- occ.dat %>% group_by(CONCEPT) %>% filter(n()>obs.n) %>% ungroup()
  f.name=paste0("all.",obs.n)
}

com <- grepl("com",type)
rar <- grepl("rar",type)

# !only carry out where we include minimum of 50 obvs! #
if (com){
  cat(type,"species selected","\n")
  
  # common rare species 
  # classify number of species at each end of the spectrum e.g. top 20 species
  n.spec= as.numeric(gsub("com.","",type))
  
  # subset occup  
  occ.cr <- data.frame(table(occ.dat$CONCEPT))
  occ.cr <- occ.cr[occ.cr$Freq>49,]
  common <- droplevels(occ.cr$Var1[order(occ.cr$Freq)][(length(occ.cr$Var1)-(n.spec-1)):length(occ.cr$Var1)])
  occ.dat <-occ.dat[occ.dat$CONCEPT%in%common,]
  f.name=type
  obs.n=49
}

# rare species
if (rar){
  cat(type,"species selected","\n")
  # common rare species 
  # classify number of species at each end of the spectrum e.g. top 20 species
  n.spec= as.numeric(gsub("rar.","",type))
  
  # subset occup  
  occ.cr <- data.frame(table(occ.dat$CONCEPT))
  occ.cr <- occ.cr[occ.cr$Freq>49,]
  rare <- droplevels(occ.cr$Var1[order(occ.cr$Freq)][1:(n.spec)])
  occ.dat <-occ.dat[occ.dat$CONCEPT%in%rare,]
  f.name=type
  obs.n=49
}

if(length(type)>1){cat(type,"species selected","\n")
                   occ.dat <-occ.dat[occ.dat$CONCEPT%in%type,]
                   f.name=paste0("ls.",length(type))}

#########################################################################################

while(!all(table(occ.dat$CONCEPT)>obs.n)){
  cat("While loop active","\n")
sp1  <-  table(occ.dat$CONCEPT)>obs.n
spp. <- names(sp1[sp1!=TRUE])
occ.dat <- occ.dat[!occ.dat$CONCEPT%in%spp.,]  

}

# this is run again as sites are removed
# derive site list  based on common sites in the covariate data and occupancy data
sites_5km <- Reduce(intersect, sites())

# subset occ.dat on common sites
occ.dat <- occ.dat[occ.dat$site_5km%in%sites_5km,]

#########################################################################
#  set up data for JAGS format

# if closure period = 2 observations are grouped together at 2-year intervals
if(closure.period==2){
  
   
   # format data
   occ <- formatOccData(taxa = occ.dat$CONCEPT,
                        site = occ.dat$TO_GRIDREF,
                        survey =  occ.dat$TO_STARTDATE,
                        closure_period = occ.dat$closure_per)
}else{
 #   use year as the closure period
   # format data
   occ <- formatOccData(taxa = occ.dat$CONCEPT,
                       site = occ.dat$TO_GRIDREF,
                       survey =  occ.dat$TO_STARTDATE)
}

# reformat to 5km grid
occ[[2]]$site_5km <- reformat_gr(occ[[2]]$site, 
                                 prec_out = 5000) # 5km2

#  create separate date frames for the occupancy info & visit info
occup <- occ[[1]]
visit <- occ[[2]]

visit <- left_join(visit[1],visit.f)


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

#########################################################################
# make sure occupancy and visit data match
if(all(occup$visit==visit$visit)){"Occupancy and visit data match"}


# use original occupancy data to compute initial values
o <- occ.dat[occ.dat$site_5km%in%visit$site_5km,]
o <- droplevels(o[o$CONCEPT%in%colnames(occup[-1]),])

# Compute observed occupancy for initial values
zobs <- table(o$CONCEPT, o$site_5km, o$closure_per)
zobs[zobs>1] <- 1
zobs <- ifelse(zobs>0,1,NA)

########################################################################
# covariates with explicit site ID need for checks
var.names = names(cov_assess)
covars_id <- list()

for (i in 1:length(var.names)){
covars_id[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=T,center=F)
names(covars_id)[i] <- var.names[[i]]
}

# double check all covs are in the right order
for (i in 1:length(covars_id)){
  if(!all(covars_id[[i]]$gr==site_id$site_name)){stop("covs are not in the right order")}
}

# covariates without site id
covars <- list()

for (i in 1:length(var.names)){
  covars[[i]] <- var_prep(cov_assess[[i]],var.x=var.x[[i]],var.names=var.names[[i]], site="gr", keep.id=F,center=T)
  names(covars)[i] <- var.names[[i]]
}
#########################################################################
# plots 
coord <- distinct(read.csv("agcensus.csv",header=T)[2:4])

id_plots <- lapply(covars_id,function(x){left_join(x,coord)})

env_plots <- lapply (id_plots,function(x1){
sy <- x1[c("E","N",start_year)]
ey <- x1[c("E","N",end_year)]

colnames(sy)[1:3] <- c("E","N","Year")
colnames(ey)[1:3] <- c("E","N","Year")

first_y <- ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data = sy, 
                               aes(x = E, y = N, fill =Year))+
  scale_fill_continuous(type = "viridis", name = start_year)

last_y <- ggplot() +
  geom_path(data = UK$britain, aes(x = long, y = lat, group = group)) + xlim(100000, 700000) +
  ylim(0, 700000)  + geom_tile(data =ey, 
                               aes(x = E, y = N, fill =Year))+
  scale_fill_continuous(type = "viridis", name = end_year)

ggpubr::ggarrange(first_y,last_y)
})

# FUNCTION ##############################################################
# check for correlation between parameters
# get combos of all variables
pair <- combn(names(covars_id),2)

# calculate correlations      
cors_mat  <- apply(pair,2,function(cname){ 
  x1 = x[[cname[1]]];y1 = x[[cname[2]]]
  
  cors <-  matrix(ncol =4,nrow=ncol(x1[-1]))
  
  for( i in 1:ncol(x1[-1])){
    
    cors[i,1] <- i
    cors[i,2]<-  round(cor(x1[i+1],y1[i+1]),2)
    cors[i,3] <-cname[1]
    cors[i,4] <-cname[2]
  }
  list(cors)
})

cors_mat
#########################################################################
# order check
if(!all(rownames(zobs[1,,])==site_id$site_name)) { cat("check sites are in the correct order")}

for (i in 1:length(covars_id)){
  if(!all(covars_id[[i]]$gr==rownames(zobs[1,,])))  {cat("check sites match covariate sites")}
}

if(!all(dimnames(zobs)[[1]]==colnames(occup[-1]))){cat("check species names are in the correct order")}

#########################################################################
# save outputs
cat("Minimum observations",obs.n,"\n")
cat(nrow(site_id),"sites","\n")
cat(ncol(occup[-1]),"species","\n")

cat(file.name <- paste0("Model_data/data_ag_",group.name,"_",f.name,"_",start_year,".",end_year,".rds"),"\n")


saveRDS(list(occup,visit,zobs,covars,closure.period,cr,cr_exc),file.name)



ggpubr::ggarrange(env_plots$RQsum,env_plots$agri,env_plots$mean_temp
                  , nrow=3)
