if (!require(CDM)) install.packages("CDM", dependencies = TRUE )
library(CDM)
library(R2jags)
library(rjags)
library(boot)
library(mcmcplots)

set.seed(12345678) #for reproducability


#####################################################Create Sampled Data#################################################################
data("data.ecpe")

# read in original data
BY04DataImport = data.ecpe$data

# create a sample from the data file:
BY04DataObs = sample(x = 1:nrow(BY04DataImport), size = nrow(BY04DataImport), replace = TRUE)

# create a new data frame using only that sample: -- Use this data frame in all analyses
BY04Data = BY04DataImport[BY04DataObs,]

BY04Data = BY04Data[,2:29]

# the q-matrix for the data is:
data.ecpe$q.matrix


################################################Begin Unidimensional 3PL Model###########################################################

nItems = ncol(BY04Data)
model01.function = function(){
  
  #measurement model specification using logistic inverse link function
  for (person in 1:N){
    for (item in 1:I){
      X[person,item] ~ dbern(c[item] + (1-c[item])*(1/(1+exp(-(a[item]*(theta[person]-b[item]))))))
    }
  }
  
 
  
  #prior distribution for theta (standardized instead of marker item)
  for (person in 1:N){
    theta[person] ~ dnorm(0,1)
  }
  
  #prior distribution for measurement model parameters
  for (item in 1:I){
    a[item] ~ dlnorm(a.mean.0,a.precision.0) #using lognormal since we are not using a marker item
    b[item] ~ dnorm(b.mean.0, b.precision.0)
    c[item] ~ dbeta(c.a.0, c.b.0)
  }
  
}

#specify hyperparameters
a.mean.0 = 0
a.variance.0 = 100
a.precision.0 = 1 / a.variance.0

b.mean.0 = 0
b.variance.0 = 100
b.precision.0 = 1 / b.variance.0

c.a.0 = 1
c.b.0 = 1 #this should give a uniform prior

#create data for JAGS to use

model01.data = list(
  N = nrow(BY04Data),
  X = BY04Data,
  I = nItems,
  a.mean.0 = a.mean.0,
  a.precision.0 = a.precision.0,
  b.mean.0 = b.mean.0,
  b.precision.0 = b.precision.0,
  c.a.0 = c.a.0,
  c.b.0 = c.b.0
)


#save model parameters
model01.parameters = c("a","b","c","theta")

#set seed for reproducable analyses
model01.seed = 15042019

#use jags.parallel() function to run simultaneous chains
model01.r2jags = jags(
  data = model01.data,
  parameters.to.save = model01.parameters,
  model.file = model01.function, 
  n.chains = 3,
  n.iter = 2000,
  n.thin = 1,
  n.burnin = 1000,
  progress.bar = "text"
)

model01.r2jags

traceplot(model01.r2jags)


###################################################Begin Multivariate 3PL Model#################################################################



# marker item:
model02.function = function(){
  # measurement model specification
  for (person in 1:N){
    X[person, 1] ~ dbern(c[1] + (1-c[1])*((1 / (1 + exp(-((a[1,1]*(theta[person,1]-b[1])) + (a[1,2]*(theta[person,2]-b[1]))))))))
    X[person, 2] ~ dbern(c[2] + (1-c[2])*((1 / (1 + exp(-(a[2,2]*(theta[person,2]-b[2])))))))
    X[person, 3] ~ dbern(c[3] + (1-c[3])*((1 / (1 + exp(-((a[3,1]*(theta[person,1]-b[3])) + (a[3,3]*(theta[person,3]-b[3]))))))))
    X[person, 4] ~ dbern(c[4] + (1-c[4])*((1 / (1 + exp(-(a[4,3]*(theta[person,3]-b[4])))))))
    X[person, 5] ~ dbern(c[5] + (1-c[5])*((1 / (1 + exp(-(a[5,3]*(theta[person,3]-b[5])))))))
    X[person, 6] ~ dbern(c[6] + (1-c[6])*((1 / (1 + exp(-(a[6,3]*(theta[person,3]-b[6])))))))
    X[person, 7] ~ dbern(c[7] + (1-c[7])*((1 / (1 + exp(-((a[7,1]*(theta[person,1]-b[7])) + (a[7,3]*(theta[person,3]-b[7]))))))))
    X[person, 8] ~ dbern(c[8] + (1-c[8])*((1 / (1 + exp(-(a[8,2]*(theta[person,2]-b[8])))))))
    X[person, 9] ~ dbern(c[9] + (1-c[9])*((1 / (1 + exp(-(a[9,3]*(theta[person,3]-b[9])))))))
    X[person,10] ~ dbern(c[10] + (1-c[10])*((1 / (1 + exp(-(a[10,1]*(theta[person,1]-b[10])))))))
    X[person,11] ~ dbern(c[11] + (1-c[11])*((1 / (1 + exp(-((a[11,1]*(theta[person,1]-b[11])) + (a[11,3]*(theta[person,3]-b[11]))))))))
    X[person,12] ~ dbern(c[12] + (1-c[12])*((1 / (1 + exp(-((a[12,1]*(theta[person,1]-b[12])) + (a[12,3]*(theta[person,3]-b[12]))))))))
    X[person,13] ~ dbern(c[13] + (1-c[13])*((1 / (1 + exp(-(a[13,1]*(theta[person,1]-b[13])))))))
    X[person,14] ~ dbern(c[14] + (1-c[14])*((1 / (1 + exp(-(a[14,1]*(theta[person,1]-b[14])))))))
    X[person,15] ~ dbern(c[15] + (1-c[15])*((1 / (1 + exp(-(a[15,3]*(theta[person,3]-b[15])))))))
    X[person,16] ~ dbern(c[16] + (1-c[16])*((1 / (1 + exp(-((a[16,1]*(theta[person,1]-b[16])) + (a[16,3]*(theta[person,3]-b[16]))))))))
    X[person,17] ~ dbern(c[17] + (1-c[17])*((1 / (1 + exp(-((a[17,2]*(theta[person,2]-b[17])) + (a[17,3]*(theta[person,3]-b[17]))))))))
    X[person,18] ~ dbern(c[18] + (1-c[18])*((1 / (1 + exp(-(a[18,3]*(theta[person,3]-b[18])))))))
    X[person,19] ~ dbern(c[19] + (1-c[19])*((1 / (1 + exp(-(a[19,3]*(theta[person,3]-b[19])))))))
    X[person,20] ~ dbern(c[20] + (1-c[20])*((1 / (1 + exp(-((a[20,1]*(theta[person,1]-b[20])) + (a[20,3]*(theta[person,3]-b[20]))))))))
    X[person,21] ~ dbern(c[21] + (1-c[21])*((1 / (1 + exp(-((a[21,1]*(theta[person,1]-b[21])) + (a[21,3]*(theta[person,3]-b[21]))))))))
    X[person,22] ~ dbern(c[22] + (1-c[22])*((1 / (1 + exp(-(a[22,3]*(theta[person,3]-b[22])))))))
    X[person,23] ~ dbern(c[23] + (1-c[23])*((1 / (1 + exp(-(a[23,2]*(theta[person,2]-b[23])))))))
    X[person,24] ~ dbern(c[24] + (1-c[24])*((1 / (1 + exp(-(a[24,2]*(theta[person,2]-b[24])))))))
    X[person,25] ~ dbern(c[25] + (1-c[25])*((1 / (1 + exp(-(a[25,1]*(theta[person,1]-b[25])))))))
    X[person,26] ~ dbern(c[26] + (1-c[26])*((1 / (1 + exp(-(a[26,3]*(theta[person,3]-b[26])))))))
    X[person,27] ~ dbern(c[27] + (1-c[27])*((1 / (1 + exp(-(a[27,1]*(theta[person,1]-b[27])))))))
    X[person,28] ~ dbern(c[28] + (1-c[28])*((1 / (1 + exp(-(a[28,3]*(theta[person,3]-b[28])))))))
  }
  
  # prior distributions for the factor: HOW DO WE "STANDARDIZE" THE FACTOR
  for (person in 1:N){
    theta[person, 1:3] ~ dmnorm(kappa[1:3], inv.phi[1:3,1:3])
  }
  
  
  # prior distribution for the factor covariance matrix
  inv.phi[1:3,1:3] ~ dwish(theta.invcov.0[1:3,1:3], theta.invcov.df.0)
  theta.cov[1:3,1:3] <- inverse(inv.phi[1:3,1:3])
  
  # fix factor means
  for (theta in 1:3){
    kappa[theta] <- 0
  }
  
  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    c[item] ~ dbeta(c.a.0, c.b.0)
    b[item] ~ dnorm(b.mean.0, b.precision.0)
  }
  
  # prior distributions for the loadings (except the first loading, which is fixed to 1.0)
  a[1,1] <- 1
  a[1,2] <- 1
  a[2,2] ~ dnorm(a.mean.0, a.precision.0)
  a[3,1] ~ dnorm(a.mean.0, a.precision.0)
  a[3,3] <- 1
  a[4,3] ~ dnorm(a.mean.0, a.precision.0) 
  a[5,3] ~ dnorm(a.mean.0, a.precision.0)
  a[6,3] ~ dnorm(a.mean.0, a.precision.0)
  a[7,1] ~ dnorm(a.mean.0, a.precision.0)
  a[7,3] ~ dnorm(a.mean.0, a.precision.0)
  a[8,2] ~ dnorm(a.mean.0, a.precision.0)
  a[9,3] ~ dnorm(a.mean.0, a.precision.0)
  a[10,1] ~ dnorm(a.mean.0, a.precision.0)
  a[11,1] ~ dnorm(a.mean.0, a.precision.0)
  a[11,3] ~ dnorm(a.mean.0, a.precision.0)
  a[12,1] ~ dnorm(a.mean.0, a.precision.0)
  a[12,3] ~ dnorm(a.mean.0, a.precision.0)
  a[13,1] ~ dnorm(a.mean.0, a.precision.0)
  a[14,1] ~ dnorm(a.mean.0, a.precision.0)
  a[15,3] ~ dnorm(a.mean.0, a.precision.0)
  a[16,1] ~ dnorm(a.mean.0, a.precision.0)
  a[16,3] ~ dnorm(a.mean.0, a.precision.0)
  a[17,2] ~ dnorm(a.mean.0, a.precision.0)
  a[17,3] ~ dnorm(a.mean.0, a.precision.0)
  a[18,3] ~ dnorm(a.mean.0, a.precision.0)
  a[19,3] ~ dnorm(a.mean.0, a.precision.0)
  a[20,1] ~ dnorm(a.mean.0, a.precision.0)
  a[20,3] ~ dnorm(a.mean.0, a.precision.0)
  a[21,1] ~ dnorm(a.mean.0, a.precision.0)
  a[21,3] ~ dnorm(a.mean.0, a.precision.0)
  a[22,3] ~ dnorm(a.mean.0, a.precision.0)
  a[23,2] ~ dnorm(a.mean.0, a.precision.0)
  a[24,2] ~ dnorm(a.mean.0, a.precision.0)
  a[25,1] ~ dnorm(a.mean.0, a.precision.0)
  a[26,3] ~ dnorm(a.mean.0, a.precision.0)
  a[27,1] ~ dnorm(a.mean.0, a.precision.0)
  a[28,3] ~ dnorm(a.mean.0, a.precision.0)
}


# specification of prior values for measurement model parameters:
#   item location
b.mean.0 = 0
b.variance.0 = 100
b.precision.0 = 1 / b.variance.0

#   item descrimination
a.mean.0 = 0
a.variance.0 = 1
a.precision.0 = 1 / a.variance.0

#values for prior c parameter
c.a.0 = 1
c.b.0 = 1

# values for prior for factor variance (based on variance of marker item) ASK ABOUT THIS
theta.cov.0 = diag(3)
theta.invcov.0 = solve(theta.cov.0)
theta.invcov.df.0 = 4 #one more than the number of factors

# next, create data for JAGS to use:
model02.data = list(
  N = nrow(BY04Data),
  X = BY04Data,
  I = ncol(BY04Data),
  b.mean.0 = b.mean.0,
  b.precision.0 = b.precision.0,
  a.mean.0 = a.mean.0,
  a.precision.0 = a.precision.0,
  c.a.0 = c.a.0,
  c.b.0 = c.b.0,
  theta.invcov.0 = theta.invcov.0,
  theta.invcov.df.0 = theta.invcov.df.0
)

model02.parameters = c("a", "b", "c", "theta.cov", "theta")

# for reproducable analyses
model02.seed = 15042019+1


#use jags.parallel() function to run simultaneous chains
model02.r2jags =  jags(
  data = model02.data,
  parameters.to.save = model02.parameters,
  model.file = model02.function, 
  n.chains = 3,
  n.iter = 5000,
  n.thin = 5,
  n.burnin = 2000,
  progress.bar = "text"
)

model02.r2jags

##############################Begin Posterior Predictive Model Fit for Multidimensional Model##################################################

# list number of simulated data sets
nSimulatedDataSets = 5000

# create one large matrix of posterior values
model02.Posterior.all = model02.r2jags$BUGSoutput$sims.matrix
dim(model02.Posterior.all)

# determine columns of posterior that go into each model matrix
aCols = grep(x = colnames(model02.Posterior.all), pattern = "a\\[")
bCols = grep(x = colnames(model02.Posterior.all), pattern = "b\\[")
cCols = grep(x = colnames(model02.Posterior.all), pattern = "c\\[")
#bText = colnames(model02.Posterior.all)[bCols]
#bCall = paste(bText, "= lambdaVec[", 1:56, "]")
#theta = matrix(data = 0, nrow = 20, ncol = 8)
covCol = grep(x = colnames(model02.Posterior.all), pattern = "theta.cov")
a = matrix(data = 0, nrow = 28, ncol=3)
b = matrix(data = 0, nrow = nrow(BY04Data), ncol = 28)
logits = matrix(data = 0, nrow = nrow(BY04Data), ncol = 28)

# save simulated covariances:
simCovModel02 = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)

# loop through data sets (can be sped up with functions and lapply)
pb = txtProgressBar()
sim = 1
for (sim in 1:nSimulatedDataSets){
  
  # draw sample from one iteration of posterior chain 
  iternum = sample(x = 1:nrow(model02.Posterior.all), size = 1, replace = TRUE)
  
  
  # get parameters for that sample: put into factor model matrices for easier generation of data
  a.vec = matrix(data = model02.Posterior.all[iternum,aCols[1:37]], ncol=1)
  b.vec = matrix(data=model02.Posterior.all[iternum,bCols], ncol=1)
  c.vec = matrix(data = model02.Posterior.all[iternum,cCols], ncol = 1)
  varTheta = matrix(data = model02.Posterior.all[iternum, covCol], nrow = 3, ncol = 3)
  
  #create proper a matrix
  a[1,1]=a.vec[1]
  a[3,1]=a.vec[2]
  a[7,1]=a.vec[3]
  a[10,1]=a.vec[4]
  a[11,1]=a.vec[5]
  a[12,1]=a.vec[6]
  a[13,1]=a.vec[7]
  a[14,1]=a.vec[8]
  a[16,1]=a.vec[9]
  a[20,1]=a.vec[10]
  a[21,1]=a.vec[11]
  a[25,1]=a.vec[12]
  a[27,1]=a.vec[13]
  a[1,2]=a.vec[14]
  a[2,2]=a.vec[15]
  a[8,2]=a.vec[16]
  a[17,2]=a.vec[17]
  a[23,2]=a.vec[18]
  a[24,2]=a.vec[19]
  a[3,3]=a.vec[20]
  a[4,3]=a.vec[21]
  a[5,3]=a.vec[22]
  a[6,3]=a.vec[23]
  a[7,3]=a.vec[24]
  a[9,3]=a.vec[25]
  a[11,3]=a.vec[26]
  a[12,3]=a.vec[27]
  a[15,3]=a.vec[28]
  a[16,3]=a.vec[29]
  a[17,3]=a.vec[30]
  a[18,3]=a.vec[31]
  a[19,3]=a.vec[32]
  a[20,3]=a.vec[33]
  a[21,3]=a.vec[34]
  a[22,3]=a.vec[35]
  a[26,3]=a.vec[36]
  a[28,3]=a.vec[37]
  
  
  # generate sample of thetas from theta distribution
  theta = rmvnorm(n = nrow(BY04Data), mean = rep(0,3), sigma = varTheta)
  
  # calculate predicted logits:
  for(i in 1:nrow(logits)){
    for(j in 1:ncol(logits)){
      logits[i,j]=a[j,1]*(theta[i,1]-b.vec[j]) + a[j,2]*(theta[i,2]-b.vec[j]) + a[j,3]*(theta[i,3]-b.vec[j])
    }
  }
  
  #simulate data
  simData = logits
  for (i in 1:ncol(logits)){
    simData[,i] = rbinom(n = nrow(logits), size = 1, prob = c.vec[i] + (1-c.vec[i]) * inv.logit(logits[,i]))
  }
  
  # calculate the value of SRMR using simulated data's covariance matrix and observed covariance matrix
  simCov = cov(simData)
  simCovModel02[sim,] = c(cov(simData))
  
  setTxtProgressBar(pb = pb, value = sim/nSimulatedDataSets)
}
close(pb)

# label values of simCor to ensure we have the right comparison
covNames = NULL
for (i in 1:ncol(simData)){
  for (j in 1:ncol(simData)){
    covNames = c(covNames, paste0("cov", i, "." , j))
  }
}
colnames(simCovModel02) = covNames

# show how one correlation compares to distribution of simulated correlations
dataCov = cov(BY04Data)
hist(simCovModel02[,1])
plot(density(simCovModel02[,1]))
lines(x = c(dataCov[1,1], dataCov[1,1]), y = c(0, max(density(simCovModel02[,1])$y)), col = 2)
quantile(simCovModel02[,1])
mean(simCovModel02[,1])
dataCov[1,1]

# create quantiles of correlations to see where each observed correlation falls
covQuantiles02 = NULL

# compute the quantiles of the observed correlations:

col = 1
for (i in 1:ncol(simData)){
  for (j in 1:ncol(simData)){
    # get empirical CDF of simulated correlation distribution
    covEcdf = ecdf(simCovModel02[,col])
    covQuantiles02 = rbind(covQuantiles02, c(i, j, summary(covEcdf), dataCov[i,j], covEcdf(dataCov[i,j])))
    
    col = col + 1
  }
}
colnames(covQuantiles02)[1:2] = c("Item 1", "Item 2")
colnames(covQuantiles02)[9:10] = c("ObsCor", "CorPctile")
covQuantiles02[which(covQuantiles02[,10] > .975 | covQuantiles02[,10] < .025),]
