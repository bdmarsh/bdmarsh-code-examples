library(rjags)
library(R2jags)
library(mvtnorm)
library(rWishart)
library(MCMCpack)



########################## Create sampled data ################################################################

#sets random number seed

set.seed(12345678) 

# read in original data
BY03DataImport = read.csv(file = "BY03Data.csv")

# create a sample from the data file:
BY03DataObs = sample(x = 1:nrow(BY03DataImport), size = nrow(BY03DataImport), replace = TRUE)

# create a new data frame using only that sample: -- Use this data frame in all analyses
BY03Data = BY03DataImport[BY03DataObs,]

# Begin Three Factor Model  ========================================================================================

# model specs:
nItems = ncol(BY03Data[paste0("Item", 1:23)])
nchains = 4
niter = 5000
nburnin = 2000
nadapt = 2000
nthin = 1

# specification of prior values for measurement model parameters:
#   item means
mu.mean.0 = 3.5 #middle of Likert scale
mu.variance.0 = 1000000 #Large to be diffuse
mu.precision.0 = 1 / mu.variance.0 #1 over variance

#   Factor loadings
lambda.mean.0 = 0
lambda.variance.0 = 1000000
lambda.precision.0 = 1 / lambda.variance.0

# unique variances
psi.df.0 = 1
psi.var.0 = apply(X = BY03Data[paste0("Item", 1:23)], MARGIN = 2, FUN = var)
psi.alpha.0 = psi.df.0 / 2
psi.beta.0 = (psi.df.0 * psi.var.0) / 2

# values for prior for factor variance (based on variance of marker item)
factor.cov.0 = diag(3)
factor.invcov.0 = solve(factor.cov.0)
factor.invcov.df.0 = 3

# factor 1: body satisfaction: 1-9
# factor 2: binge eating: 10-16
# factor 3: thinness: 17-23
model01a.function = function(){
  
  # measurement model specification
  for (person in 1:N){
    mean[person,  1] <- mu[1]  + lambda[1,  1]*xfactor[person, 1]
    mean[person,  2] <- mu[2]  + lambda[2,  1]*xfactor[person, 1]
    mean[person,  3] <- mu[3]  + lambda[3,  1]*xfactor[person, 1]
    mean[person,  4] <- mu[4]  + lambda[4,  1]*xfactor[person, 1]
    mean[person,  5] <- mu[5]  + lambda[5,  1]*xfactor[person, 1]
    mean[person,  6] <- mu[6]  + lambda[6,  1]*xfactor[person, 1]
    mean[person,  7] <- mu[7]  + lambda[7,  1]*xfactor[person, 1]
    mean[person,  8] <- mu[8]  + lambda[8,  1]*xfactor[person, 1]
    mean[person,  9] <- mu[9]  + lambda[9,  1]*xfactor[person, 1]
    mean[person, 10] <- mu[10] + lambda[10, 2]*xfactor[person, 2]
    mean[person, 11] <- mu[11] + lambda[11, 2]*xfactor[person, 2]
    mean[person, 12] <- mu[12] + lambda[12, 2]*xfactor[person, 2]
    mean[person, 13] <- mu[13] + lambda[13, 2]*xfactor[person, 2]
    mean[person, 14] <- mu[14] + lambda[14, 2]*xfactor[person, 2]
    mean[person, 15] <- mu[15] + lambda[15, 2]*xfactor[person, 2]
    mean[person, 16] <- mu[16] + lambda[16, 2]*xfactor[person, 2]
    mean[person, 17] <- mu[17] + lambda[17, 3]*xfactor[person, 3]
    mean[person, 18] <- mu[18] + lambda[18, 3]*xfactor[person, 3]
    mean[person, 19] <- mu[19] + lambda[19, 3]*xfactor[person, 3]
    mean[person, 20] <- mu[20] + lambda[20, 3]*xfactor[person, 3]
    mean[person, 21] <- mu[21] + lambda[21, 3]*xfactor[person, 3]
    mean[person, 22] <- mu[22] + lambda[22, 3]*xfactor[person, 3]
    mean[person, 23] <- mu[23] + lambda[23, 3]*xfactor[person, 3]
    
    for (item in 1:I){
      X[person, item] ~ dnorm(mean[person,item], inv.psi[item])  
    }
  }
  
  # prior distributions for the factor:
  for (person in 1:N){
    xfactor[person, 1:3] ~ dmnorm(kappa, inv.phi[,])
  }
  
  # prior distribution for the factor covariance matrix
  inv.phi ~ dwish(factor.invcov.0, factor.invcov.df.0)
  
  # fix factor means
  for (factor in 1:3){
    kappa[factor] <- 0
  }

  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
    inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0[item])
  }
  
  # prior distributions for the loadings (except the first loading of each factor, which is fixed to 1.0)
  lambda[1,1] <- 1
  lambda[2,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[3,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[4,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[5,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[6,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[7,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[8,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[9,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[10,2] <- 1
  lambda[11,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[12,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[13,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[14,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[15,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[16,2] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[17,3] <- 1
  lambda[18,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[19,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[20,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[21,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[22,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[23,3] ~ dnorm(lambda.mean.0, lambda.precision.0)
  
  # saved parameters
  for (item in 1:I){
    psi[item] <- 1/inv.psi[item]
  }
  phi = inverse(inv.phi)
  
}

# for reproducable seeds (without parallel JAGS)
model01a.seed = 12345678

#   Note: here, the random number seed cannot be the same per seed or the chains will be the same
RNGname = c("Wichmann-Hill","Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")
if (RNGname[1] %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                      "Super-Duper", "Mersenne-Twister")) {
  RNGname[1] <- paste("base::", RNGname[1], sep = "")
}

model01a.init.values <- vector("list", nchains)
for (i in 1:nchains) {
  model01a.init.values[[i]]$.RNG.name <- RNGname[1]
  model01a.init.values[[i]]$.RNG.seed <- model01a.seed + i
}
  
# next, create data for JAGS to use:
model01a.data = list(
  N = nrow(BY03Data),
  X = BY03Data[paste0("Item", 1:23)],
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  psi.alpha.0 = psi.alpha.0,
  psi.beta.0 = psi.beta.0,
  factor.invcov.0 = factor.invcov.0,
  factor.invcov.df.0 = factor.invcov.df.0
)

#model parameters for jags
model01a.parameters = c("mu", "lambda",  "psi", "phi", "xfactor")


# run model
model01a.r2jags =  jags(
  data = model01a.data,
  inits = model01a.init.values,
  parameters.to.save = model01a.parameters,
  model.file = model01a.function,
  n.chains = nchains,
  n.iter = niter,
  n.burnin = nburnin,
  n.thin = nthin
)

#display model
model01a.r2jags

# examining the factor scores
dim(model01a.r2jags$BUGSoutput$sims.matrix)
colnames(model01a.r2jags$BUGSoutput$sims.matrix)

# Begin Posterior Predictive Model Checks ==============================================================================

# list number of simulated data sets
nSimulatedDataSets = 5000

# create one large matrix of posterior value by disentangling chains
model01a.Posterior.all = model01a.r2jags$BUGSoutput$sims.matrix


# determine columns of posterior that go into each model matrix
muCols =grep(x = colnames(model01a.Posterior.all), pattern = "mu")
colnames(model01a.r2jags$BUGSoutput$sims.matrix)[muCols]
lambdaCols = grep(x = colnames(model01a.Posterior.all), pattern = "lambda")
colnames(model01a.r2jags$BUGSoutput$sims.matrix)[lambdaCols]
psiCols = grep(x = colnames(model01a.Posterior.all), pattern = "psi")
colnames(model01a.r2jags$BUGSoutput$sims.matrix)[psiCols]
phiCols = grep(x = colnames(model01a.Posterior.all), pattern = "phi")
colnames(model01a.r2jags$BUGSoutput$sims.matrix)[phiCols]

# save simulated correlations:
simCor = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*(nItems-1)/2)
simCovModel01a = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)
simSRMR = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)
simRMSEA = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)

# save model-based covariances:
model01aCov = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)

# model DF
modelDF = (nItems*(nItems+1)/2) - 46

dataCov = cov(BY03Data[paste0("Item", 1:23)])
detDataCov = det(dataCov)

# loop through data sets (can be sped up with functions and lapply)
pb = txtProgressBar()
sim = 1
for (sim in 1:nSimulatedDataSets){
  
  # draw sample from one iteration of posterior chain 
  iternum = sample(x = 1:nrow(model01a.Posterior.all), size = 1, replace = TRUE)
  
  # get parameters for that sample: put into factor model matrices for easier generation of data
  mu = matrix(data = model01a.Posterior.all[iternum, muCols], ncol = 1)
  psi = diag(model01a.Posterior.all[iternum, psiCols])
  phi = matrix(model01a.Posterior.all[iternum, phiCols], nrow = 3, ncol = 3)
  
  lambda = matrix(data = 0, ncol = 3, nrow = 23)
  lambda[1,1] = model01a.Posterior.all[iternum, lambdaCols[1]]
  lambda[2,1] = model01a.Posterior.all[iternum, lambdaCols[2]]
  lambda[3,1] = model01a.Posterior.all[iternum, lambdaCols[3]]
  lambda[4,1] = model01a.Posterior.all[iternum, lambdaCols[4]]
  lambda[5,1] = model01a.Posterior.all[iternum, lambdaCols[5]]
  lambda[6,1] = model01a.Posterior.all[iternum, lambdaCols[6]]
  lambda[7,1] = model01a.Posterior.all[iternum, lambdaCols[7]]
  lambda[8,1] = model01a.Posterior.all[iternum, lambdaCols[8]]
  lambda[9,1] = model01a.Posterior.all[iternum, lambdaCols[9]]
  lambda[10,2] = model01a.Posterior.all[iternum, lambdaCols[10]]
  lambda[11,2] = model01a.Posterior.all[iternum, lambdaCols[11]]
  lambda[12,2] = model01a.Posterior.all[iternum, lambdaCols[12]]
  lambda[13,2] = model01a.Posterior.all[iternum, lambdaCols[13]]
  lambda[14,2] = model01a.Posterior.all[iternum, lambdaCols[14]]
  lambda[15,2] = model01a.Posterior.all[iternum, lambdaCols[15]]
  lambda[16,2] = model01a.Posterior.all[iternum, lambdaCols[16]]
  lambda[17,3] = model01a.Posterior.all[iternum, lambdaCols[17]]
  lambda[18,3] = model01a.Posterior.all[iternum, lambdaCols[18]]
  lambda[19,3] = model01a.Posterior.all[iternum, lambdaCols[19]]
  lambda[20,3] = model01a.Posterior.all[iternum, lambdaCols[20]]
  lambda[21,3] = model01a.Posterior.all[iternum, lambdaCols[21]]
  lambda[22,3] = model01a.Posterior.all[iternum, lambdaCols[22]]
  lambda[23,3] = model01a.Posterior.all[iternum, lambdaCols[23]]
  
  
  # create model-implied mean and covariance matrix (marginal for X)
  meanVec = mu
  covMat = lambda %*% phi %*% t(lambda) + psi
  model01aCov[sim, ] = c(covMat)
  
  # randomly draw data with same sample size from MVN with mean=meanVec and cov=covMat
  simData = rmvnorm(n = nrow(BY03Data), mean = meanVec, sigma = covMat)
  
  # create sample statistics from simulated data (we'll use correlation matrix, starting with upper triangle)
  simCor[sim,] = matrix(data = c(cor(simData)[upper.tri(cor(simData))]), nrow = 1)
  
  # calculate the value of SRMR using simulated data's covariance matrix and observed covariance matrix
  simCov = cov(simData)
  simCovModel01a[sim,] = c(cov(simData))
  difCov = dataCov-simCov
  stdDifCov = sqrt(solve(diag(diag(dataCov)))) %*% difCov %*% sqrt(solve(diag(diag(dataCov))))
  
  # using formula from book:
  simSRMR[sim,1] = sum(colSums(stdDifCov*stdDifCov))/((ncol(simCov)*(ncol(simCov)-1))/2)
  
  # can also do a similar process for RMSEA (assuming covariance matrix is from ML estimation - discrepancy function)
  simRMSEA[sim,1] = sqrt((log(det(simCov)) - log(det(dataCov)) + sum(diag(dataCov %*% solve(simCov))) - nItems)/modelDF)
  
  setTxtProgressBar(pb = pb, value = sim/nSimulatedDataSets)
}
close(pb)

# first, we examine the posterior predictive distribution of SRMR (p. 241)
hist(simSRMR[,1], xlab = "Simulated SRMR", main = "SRMR Histogram")
plot(density(simSRMR))
quantile(simSRMR)
mean(simSRMR)

# next we can examine the posterior predictive distribution of RMSEA
hist(simRMSEA[,1], xlab = "Simulated RMSEA", main = "RMSEA Histogram")
plot(density(simRMSEA))
quantile(simRMSEA)
mean(simRMSEA)

# label values of simCor to ensure we have the right comparison
corNames = NULL
for (i in 1:(ncol(simData)-1)){
  for (j in (i+1):ncol(simData)){
    corNames = c(corNames, paste0("cor", i, "." , j))
  }
}
colnames(simCor) = corNames


# show how one correlation compares to distribution of simulated correlations
dataCor = cor(BY03Data[paste0("Item", 1:23)])
hist(simCor[,1])
plot(density(simCor[,1]))
lines(x = c(dataCor[1,2], dataCor[1,2]), y = c(0, 5), col = 2)
quantile(simCor[,1])
mean(simCor[,1])

# create quantiles of correlations to see where each observed correlation falls
corQuantiles = NULL

# compute the quantiles of the observed correlations:

col = 1
for (i in 1:(ncol(simData)-1)){
  for (j in (i+1):ncol(simData)){
    # get empirical CDF of simulated correlation distribution
    corEcdf = ecdf(simCor[,col])
    corQuantiles = rbind(corQuantiles, c(i, j, summary(corEcdf), dataCor[i,j], corEcdf(dataCor[i,j])))
    
    col = col + 1
  }
}
colnames(corQuantiles)[1:2] = c("Item 1", "Item 2")
colnames(corQuantiles)[9:10] = c("ObsCor", "CorPctile")
corQuantiles[which(corQuantiles[,10] > .975 | corQuantiles[,10] < .025),]

# Model 2a: Saturated Mean/Variance/Covariance Model ==================================================================

# model specs:
nItems = ncol(BY03Data[paste0("Item", 1:23)])
nchains = 4
niter = 5000
nburnin = 2000
nadapt = 2000
nthin = 1


# prior values for each mean -- joinly MVN (but diagonal covariance/precision matrix)
mean.0 = 3.5
mean.variance.0 = 1000000
mean.precision.0 = 1/mean.variance.0
meanVec.0 = matrix(data = rep(mean.0, nItems), nrow = nItems, ncol = 1)
meanPrecision.0 = diag(nItems)*mean.precision.0

# prior values for variances -- variances of observed data
R.0 = (apply(X = BY03Data[paste0("Item", 1:23)], MARGIN = 2, FUN = var) * diag(nItems))/nItems
k.0 = nItems


# load data into list
model02a.data = list(N = nrow(BY03Data),
                     X = BY03Data[paste0("Item", 1:23)],
                     I = nItems, 
                     meanVec.0 = meanVec.0,
                     meanPrecision.0 = meanPrecision.0,
                     R.0 = R.0,
                     k.0 = k.0)

# save parameters
model02a.parameters = c("mu", "sigma", "deviance")

# copy initial values
model02a.init.values = model01a.init.values

#R2Jags model

model02a.function = function(){
  # define model likelihood as multivariate normal
  for (person in 1:N){
    X[person, 1:I] ~ dmnorm(mu[1:I], sigmainv[1:I, 1:I])
  }
  
  # prior for mean vector
  mu[1:I] ~ dmnorm(meanVec.0, meanPrecision.0) 
  
  # prior for inverse covariance matrix
  sigmainv[1:I, 1:I] ~ dwish(R.0, k.0)
  
  # save covariance matrix 
  sigma = inverse(sigmainv)
}


model02a.r2jags =  jags(
  data = model02a.data,
  inits = model02a.init.values,
  parameters.to.save = model02a.parameters,
  model.file = model02a.function,
  n.chains = nchains,
  n.iter = niter,
  n.burnin = nburnin,
  n.thin = nthin
)
model02a.r2jags

#########################################################


#view traceplots
traceplot(model01a.r2jags)

#view histograms
for (i in 1:23) {
  hist(BY03Data[,i], xlab= "Item Score", main=paste0("Histogram of Item ",i))
  
}


