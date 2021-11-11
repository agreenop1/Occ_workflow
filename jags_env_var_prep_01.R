library(dplyr)
# FUNCTION ##########################################################################
# subset yearly variable data either by mean (smry=M), first year of closure period (F) or
# last year of closure period (L)
var_subset <- function(x,var.x,smry){
  
  # take mean of two years
  if(smry=="M"){
    
    s <- sym(var.x) # cannot pass quote to dplyr
    s=enquo(s)
    x.  <- x %>% group_by(gr,closure_per) %>% summarise(mean.x= mean(!!s)) # mean between years
    x1  <- merge(x.,clsr_per[clsr_per$YEAR%in%seq(start_year,end_year,2),])
    
    # first year of closure period used
  }else if(smry=="F") {
    
    x. <- x[x$YEAR%in%seq(start_year,end_year,2),]  
    x1 <- data.frame(x.["closure_per"],x.["gr"], x.[var.x],x.["YEAR"])
    
    #second year of closure period used
  }else {
    
    x. <- x[x$YEAR%in%seq(start_year+1,end_year+1,2),]  
    x1 <- data.frame(x.["closure_per"],x.["gr"], x.[var.x],x.["YEAR"])
    
  }
  return(x1)
}

env_cov <- read.csv("covars.csv")[-1] # environmental data
env_cov$YT_ID <- NULL

# closure period
closure.period=2
lag = T

# years
#  if function to ensure the correct occurrence records are selected across years
start_year =1994; end_year = 2010

years <- if(lag&closure.period==1){start_year:(closure.period+end_year)
}else if(lag&closure.period==2) {start_year:(closure.period+end_year+1)
}else if(!lag&closure.period==1) {(start_year-closure.period):(end_year)
}else {(start_year-closure.period):(end_year+1)}


# calc closure
n_periods <- rep(1:(length(years)/closure.period),each=2) # Create vector half the length of years
clsr_per <- data.frame(n_periods,years) #  assign new temporal ID to each year
colnames(clsr_per)[1:2] <- c("closure_per","YEAR")


##############################################
# subset environmental covars on sites & years
env_cov <- read.csv("covars.csv")[-1] # environmental data
env_cov$YT_ID <- NULL
env <- env_cov
colnames(env)[1:2] <- c("gr", "YEAR")

# summary env
summary.env="M" # mean

# merge with new closure periods
env <- left_join(env,clsr_per)

# add all semi nat data
env$semi <- env$Grassland + env$Tree_shrub +env$vegetation_brackish_water

# temperature
temp <- var_subset(env,"temp_anomaly",smry=summary.env)
colnames(temp)[4] <- "year" 

# agriculture 
agri <- var_subset(env,"Agriculture",smry=summary.env)
colnames(agri)[4] <- "year" 

# all semi
semi <- var_subset(env,"semi",smry=summary.env)
colnames(semi)[4] <- "year"

env_covar <- list(temp=temp,agri=agri,semi=semi)
saveRDS(env_covar,"env.cov.list.rds")
