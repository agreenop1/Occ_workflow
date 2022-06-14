# Practical 1 - building a model from BUGS/JAGS code

library(nimble)


##########################
### Simulate some data ###
##########################

# Choose sample sizes and prepare observed data array y
set.seed(24) # So we all get same data set
M <- 10 # Number of sites
J <- 2 # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data
# Parameter values
psi <- 0.8 # Probability of occupancy or presence
p <- 0.5 # Probability of detection
# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli
# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

######################################################
##### Define and create an operational model ########
######################################################

# Build a model in R directly using nimbleCode

occNimble <- nimbleCode({
  for (i in 1:nsite) {
    z[i] ~ dbern(psi) # True occupancy status
    for (j in 1:nvisit) {
      y[i,j] ~ dbern(z[i] * p) # Observed data
    }
  }
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  sum.z <- sum(z[1:nsite])
})


# Create a model object (uncompiled)
# data and initial values can be declared here or later on using the setData() and setInits()

occModel <- nimbleModel(occNimble, constants = list(nsite = M, nvisit = J))


# read a pre existing one and creates a model object

occJAGS<- readBUGSmodel('Occ_model_simple', data =  list(y = y, M = nrow(y), J = ncol(y)),
                        inits = list(z = apply(y, 1, max)))



##########################################
##### Manipulate operational models ######
##########################################

# NIMBLE models are objects you can query and manipulate, you can use it to set or get values


# Determining the nodes and variables in a model

occJAGS$getVarNames()

#note that lifted nodes can also correspond to different parametrizations

occJAGS$getNodeNames()


# Model variables can be accessed and set just as in R using $ and [[ ]].

occJAGS$psi<-0.8
occJAGS[["psi"]]


occJAGS$y[1:4,1:2]
head(occJAGS$y)
# Each node can be accessed it directly

occJAGS$z[1]
occJAGS[["z[1]"]]<-0
occJAGS[["z[1:3]"]]

# Using node names enables functions in R or NIMBLE that access any single node by a name to be written,
# regardless of the dimensionality of the variable in which it is embedded



# Checking if a node holds data | you can query whether a node is flagged as data


occJAGS$isData('z[1]')
occJAGS$isData('z')
occJAGS$isData('y')

# now for the nimble-written model

occModel$isData('y')


##############################################
#### setting new data and initial values #####
##############################################

# set data and initials for the uncompiled model 
# usually you will do this with the compile version or directly in the nimbleModel function

occModel$setData(list(y = y))
occModel$setInits(list(z = apply(y, 1, max),psi = 0.7, p = 0.33))
occModel$isData('y')

# Query the model's relationships

occModel$getDependencies(c('psi'))

library(igraph)

plot(occModel$getGraph())


#############################################
##### Simulating  directly from a model #####
#############################################


# NIMBLE models can be used to simulate data 
# Note that constants and initial values for parameters are needed for simulation


occModelsim <- nimbleModel(occNimble, constants = list(nsite = 100, nvisit = 4),
                           inits = list(psi = 0.7, p = 0.33))

# y and z nodes are not initialized

# If we try to calculate the log probabilities of these nodes at this point, we'll get NA.

occModelsim$calculate('y')

# simulate() to specify which nodes to simulate e.g. 

occModelsim$simulate("y")

# error are because y depends on z which is not initialized.
# we simulate z before y to get aorund this issue
# Use getDependencies to get all nodes needing to be simulated given the specified parameters

nodesToSim <- occModelsim$getDependencies(c("psi", "p"),
                                       self = F, downstream = T)

nodesToSim
#  self = FALSE removes psi or p in the return vector. 
#  downstream = TRUE recursively looks at the dependencies of nodes instead of
#  the nodes that depend directly on psi and p


occModelsim$simulate(nodesToSim)
occModelsim$y
occModelsim$z

occModelsim$psi
occModelsim$p


# It is better to work with the compiled model!


CoccModelsim <- compileNimble(occModelsim) # note that we compile the operational model not the code!
nodesToSim <- CoccModelsim$getDependencies(c("psi", "p"),
                                        self = F, downstream = T)
CoccModelsim$simulate(nodesToSim)

CoccModelsim$y


# NIMBLE's simulate functions  do not overwrite data values by default.
# Values of data variables can be set using setData() or replaced using resetData() to inform NIMBLE which
# nodes should be treated as data

CoccModelsim$isData('y')

# retrieve simulated data
y.sim <- CoccModelsim$y


CoccModelsim$setData(list(y = y.sim))
CoccModelsim$isData('y')

#however if we had some missing values or a mix of data/non data we would have to use the reset function
# to declare which nodes should be flagged as data.


### Missing values ####

# simulate some data

set.seed(123)



# we can make a copy of the model

simOcc_miss <- CoccModelsim$newModel()
simOcc_miss$resetData()

y.sim[3,1:2]<-NA
simOcc_miss$setData(list(y=y.sim))
simOcc_miss$isData('y')


configureMCMC(simOcc_miss)


##########################
### End of practical 1 ###
##########################
