
set.seed(12345678)


library(AICcmodavg)
library(rjags)
library(R2jags)
library(CDM)
library(ggplot2)
library(mcmc)
library(mcmcplots)
library(boot)

data("fraction.subtraction.data")

# read in original data
BY05DataImport1 = fraction.subtraction.data

# create a sample from the data file:
BY05DataObs2 = sample(x = 1:nrow(BY05DataImport1), size = nrow(BY05DataImport1), replace = TRUE)

# create a new data frame using only that sample: -- Use this data frame in all analyses
BY05Data1 = BY05DataImport1[BY05DataObs2,]

#specify number of items
nItems1 = ncol(BY05Data1)

#######################################Begin 3PNO Model##########################################################

#measurment model specifications
model3pl.function = function(){
  for (person in 1:N){
    for (item in 1:I){
      X[person, item] ~ dbern(c[item] + (1-c[item])*phi(mu[item] +lambda[item]*theta[person]))
    }
  }
  
  #prior distribution for ability
  for(person in 1:N){
    theta[person] ~ dnorm(0,1)
  }
  
  #prior distribution item slope and intercept
  for(item in 1:I){
    mu[item] ~ dnorm(0,1)
    lambda[item] ~ dlnorm(0,2)
    c[item] ~ dbeta(1,1)
  }
}

model3pl.data = list(
  N = nrow(BY05Data1),
  X = BY05Data1,
  I = nItems1
)

model3pl.parameters = c("mu", "lambda", "theta", "c")

model3pl.r2jags = jags(
  data = model3pl.data,
  parameters.to.save = model3pl.parameters,
  model.file = model3pl.function,
  n.chains = 3,
  n.iter = 5000,
  n.thin = 1,
  n.burnin = 2000,
  progress.bar = "text"
)

model3pl.r2jags

DIC(model3pl.r2jags)

#####################################Begin 3P-RH Model##############################################################

#measurment model specifications
model3prh.function = function(){
  for (person in 1:N){
    for (item in 1:I){
      X[person, item] ~ dbern(phi((sqrt(1 + exp(-delta[item]*theta[person]))*(lambda[item]*theta[person] + mu[item])) / sqrt(2)))
    }
  }
  
  #prior distribution for ability
  for(person in 1:N){
    theta[person] ~ dnorm(0,1)
  }
  
  #prior distribution item slope and intercept
  for(item in 1:I){
    mu[item] ~ dnorm(0,1)
    lambda[item] ~ dlnorm(0,2)
    delta[item] ~ dnorm(0,1)
  }
}

model3prh.data = list(
  N = nrow(BY05Data1),
  X = BY05Data1,
  I = nItems1
)

model3prh.parameters = c("mu", "lambda", "theta", "delta")

model3prh.r2jags = jags(
  data = model3prh.data,
  parameters.to.save = model3prh.parameters,
  model.file = model3prh.function,
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 2000,
  progress.bar = "text"
)

model3prh.r2jags

####################################### Begin Part 2 ###########################################

# read in original data
BY05DataImport = read.csv("conspiracies.csv")

# create a sample from the data file:
BY05DataObs = sample(x = 1:nrow(BY05DataImport), size = nrow(BY05DataImport), replace = TRUE)

# create a new data frame using only that sample: -- Use this data frame in all analyses
BY05Data = BY05DataImport[BY05DataObs,]

BY05Data = BY05Data[,1:10]

# subtract 1 from every item, making them range from 0-4
for (i in 1:10){
  BY05Data[,i] = BY05Data[,i]-1
}



#################################Begin Binomial Model####################################################


nItems = ncol(BY05Data)
model01.function = function(){
  
  #measurement model specification using logistic inverse link function
  for (person in 1:N){
    for (item in 1:I){
      X[person,item] ~ dbin(max(min(exp(-(a[item]*(theta[person]-b[item]))),0.9999),0.0001),4)
    }
  }
  
  #prior distribution for theta
  for (person in 1:N){
    theta[person] ~ dnorm(0,1)
  }
  
  
  #prior distribution for discrimination parameter
  for (item in 1:I){
    a[item] ~ dnorm(a.mean.0,a.precision.0);T(0,)
  }
  
  for(item in 1:I){
    b[item] ~ dnorm(b.mean.0,b.precision.0)
  }
  
  
}

# specification of prior values for measurement model parameters:

#item discriminations
a.mean.0 = 0
a.variance.0 = 100
a.precision.0 = 1 / a.variance.0

#item locations
b.mean.0 = 0
b.variance.0 = 100
b.precision.0 = 1 / b.variance.0


#create data for JAGS to use

model01.data = list(
  N = nrow(BY05Data),
  X = BY05Data,
  I = nItems,
  a.mean.0 = a.mean.0,
  a.precision.0 = a.precision.0,
  b.mean.0 = b.mean.0,
  b.precision.0 = b.precision.0
)

#save model parameters
model01.parameters = c("a","b","theta")

#set seed for reproducable analyses
model01.seed = 15042019


#use jags.parallel() function to run simultaneous chains
model01.r2jags = jags(
  data = model01.data,
  parameters.to.save = model01.parameters,
  model.file = model01.function, 
  n.chains = 3,
  n.iter = 3000,
  n.thin = 1,
  n.burnin = 2000,
  progress.bar = "text"
)

model01.r2jags

#################################### Begin Graded Response Model#####################################


#re-transform the data
BY05DataGRM = BY05Data

for (i in 1:10){
  BY05DataGRM[,i] = BY05Data[,i]+1
}

# marker item:
modelgr.function = function(){
  
  # measurement model specification
  for (person in 1:N){
    for (item in 1:I){
      
      # form cumulative probability item response functions
      CProb[person, item, 1] <- 1
      for (cat in 2:C[item]){
        CProb[person, item, cat] <- phi(a[item]*(theta[person]-b[item, (cat-1)]))  
      }
      
      # form probability response is equal to each category
      for (cat in 1:(C[item] - 1)){
        Prob[person, item, cat] <- CProb[person, item, cat] - CProb[person, item, cat+1]
      }
      Prob[person, item, C[item]] <- CProb[person, item, C[item]]
      
      X[person, item] ~ dcat(Prob[person, item, 1:C[item]])
    }
  }
  
  # prior distributions for the factor:
  for (person in 1:N){
    theta[person] ~ dnorm(0, 1)
  }
  
  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    
    # create parameters that are unbounded, then sort
    for (cat in 1:(C[item]-1)){
      b.star[item, cat] ~ dnorm(b.mean.0, b.precision.0)  
    }
    b[item, 1:(C[item]-1)] <- sort(b.star[item, 1:(C[item]-1)])
    
    # loadings are set to be all positive
    a[item] ~ dnorm(a.mean.0, a.precision.0);T(0,)
    
  }
  
}

# specification of prior values for measurement model parameters:
#   item intercepts
b.mean.0 = 0
b.variance.0 = 100
b.precision.0 = 1 / b.variance.0

#   Factor loadings -- these are the discriminations
a.mean.0 = 0
a.variance.0 = 100
a.precision.0 = 1 / a.variance.0

# next, create data for JAGS to use:
modelgr.data = list(
  N = nrow(BY05DataGRM),
  X = BY05DataGRM,
  C = unlist(apply(X = BY05DataGRM[,1:10], MARGIN = 2, FUN = max)),
  I = 10,
  b.mean.0 = b.mean.0,
  b.precision.0 = b.precision.0,
  a.mean.0 = a.mean.0,
  a.precision.0 = a.precision.0
)

modelgr.init = function(){
  list("a" = runif(10, 1, 2),
       "b.star" = cbind(rep(1, 10), rep(0, 10), rep(-1, 10), rep(-2, 10)))
}

modelgr.parameters = c("a", "b",  "theta")



modelgr.r2jags =  jags(
  data = modelgr.data,
  inits = modelgr.init,
  parameters.to.save = modelgr.parameters,
  model.file = modelgr.function,
  n.chains = 3,
  n.iter = 3000,
  n.thin = 1,
  n.burnin = 2000,
  progress.bar = "text"
)

modelgr.r2jags

##################################### Begin Generalized Partial Credit Model ########################################


# marker item:
modelgpc.function = function(){
  
  # measurement model specification
  for (person in 1:N){
    for (item in 1:I){
      
      for (cat in 1:C[I]){
        eta[person, item, cat] <- a[item] * (theta[person] - b[item, cat])
        psum[person, item, cat] <- sum(eta[person, item, 1:cat])
        exp.psum[person, item, cat] <- exp(psum[person, item, cat])
        prob[person, item, cat] <- exp.psum[person, item, cat]/sum(exp.psum[person, item, 1:C[item]])
      }
      
      X[person, item] ~ dcat(prob[person, item, 1:C[item]])
    }
  }
  
  # prior distributions for the factor:
  for (person in 1:N){
    theta[person] ~ dnorm(0, 1)
  }
  
  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    
    b[item, 1] <- 0
    
    # create parameters that are unbounded, then sort
    for (cat in 2:C[item]){
      b[item, cat] ~ dnorm(b.mean.0, b.precision.0)  
    }
    
    # loadings are set to be all positive
    a[item] ~ dnorm(a.mean.0, a.precision.0);T(0,)
    
  }
  
}


# specification of prior values for measurement model parameters:
#   item intercepts
b.mean.0 = 0
b.variance.0 = 100
b.precision.0 = 1 / b.variance.0

#   Factor loadings -- these are the discriminations
a.mean.0 = 0
a.variance.0 = 100
a.precision.0 = 1 / a.variance.0

# next, create data for JAGS to use:
modelgpc.data = list(
  N = nrow(BY05DataGRM),
  X = BY05DataGRM,
  C = unlist(apply(X = BY05DataGRM[,1:10], MARGIN = 2, FUN = max)),
  I = 10,
  b.mean.0 = b.mean.0,
  b.precision.0 = b.precision.0,
  a.mean.0 = a.mean.0,
  a.precision.0 = a.precision.0
)

modelgpc.init = function(){
  list("a" = runif(10, 1, 2),
       "b" = cbind(rep(NA, 10), rep(1, 10), rep(0, 10), rep(-1, 10), rep(-2, 10)))
}

modelgpc.parameters = c("a", "b","theta")


modelgpc.r2jags =  jags(
  data = modelgpc.data,
  inits = modelgpc.init,
  parameters.to.save = modelgpc.parameters,
  model.file = modelgpc.function,
  n.chains = 3,
  n.iter = 3000,
  n.thin = 1,
  n.burnin = 1000,
  progress.bar = "text"
)

modelgpc.r2jags


