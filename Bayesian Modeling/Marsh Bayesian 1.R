
#loads necessary packages
library(rjags)
library(HDInterval)

set.seed(12345678)

# read in original data
HW02DataImport = read.csv(file = "BY02Data.csv")

# create a sample from the data file:
BY02DataObs = sample(x = 1:nrow(BY02DataImport), size = nrow(BY02DataImport), replace = TRUE)

# create a new data frame using only that sample: -- Use this data frame in all analyses
BY02Data = BY02DataImport[BY02DataObs,]



########Begin fit of first model########

#specifies the linear model
FullModel = lm(score~exp4+enthusiasm+CvMG+CvRG+CvBG+exp4:CvMG+exp4:CvRG+exp4:CvBG+enthusiasm:CvMG+enthusiasm:CvRG+enthusiasm:CvBG
               +exp4:enthusiasm+exp4:enthusiasm:CvMG+exp4:enthusiasm:CvRG+exp4:enthusiasm:CvBG, data = BY02Data)

# specify specs for analysis
n.chains = 4

# specify values of hyperparameters for prior distributions


# for betas:
betaMean.0 = 0
betaVariance.0 = 100000000
betaTau.0 = 1/betaVariance.0

# take data frame and create model matrix -- which gives us the columns of X:
BayesianFullModel.Matrix = model.matrix(object = FullModel)

# append with weight and then select only columns needed for model:
BayesianFullModel.Matrix = data.frame(cbind(BY02Data$score, BayesianFullModel.Matrix))
names(BayesianFullModel.Matrix) = c("score", "intercept", "exp4", "enthusiasm", "CvMG", "CvRG", "CvBG", "exp4*CvMG", "exp4*CvRG", "exp4*CvBG", "enthusiasm*CvMG", "enthusiasm*CvRG", "enthusiasm*CvBG", "exp4*enthusiasm", "exp4*enthusiasm*CvMG", "exp4*enthusiasm*CvRG", "exp4*enthusiasm*CvBG")
BayesianFullModel.Matrix = BayesianFullModel.Matrix[c("score", "exp4", "enthusiasm", "CvMG", "CvRG", "CvBG")]

# add data and N to list object that we will pass to JAGS
BayesianFullModel.JAGS.Data = list(N = nrow(BY02Data),
                                   score = BayesianFullModel.Matrix$score,
                                   exp4 = BayesianFullModel.Matrix$exp4,
                                   enthusiasm = BayesianFullModel.Matrix$enthusiasm,
                                   CvMG = BayesianFullModel.Matrix$CvMG,
                                   CvRG = BayesianFullModel.Matrix$CvRG,
                                   CvBG = BayesianFullModel.Matrix$CvBG,
                                   betaMean.0 = betaMean.0,
                                   betaTau.0 = betaTau.0)


BayesianFullModel.Syntax = "

model{


# prior distributions: 
#   Note that terms on the left are model parameter names we create here
beta.0             ~ dnorm(betaMean.0, betaTau.0) # prior for the intercept
beta.exp4          ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional slope of exp4
beta.enth          ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional slope of enthusiasm
beta.MG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between MG group and Control
beta.RG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between RG group and Control
beta.BG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between BG group and Control
beta.exp4.MG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and MG group
beta.exp4.RG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and RG group
beta.exp4.BG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and BG group
beta.enth.MG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and MG group
beta.enth.RG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and RG group
beta.enth.BG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and BG group
beta.exp4.enth     ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and enthusiasm
beta.exp4.enth.MG  ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4, enthusiasm, and MG group
beta.exp4.enth.RG  ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4, enthusiasm, and RG group
beta.exp4.enth.BG  ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4, enthusiasm, and BG group
gamma.0            ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of residual precision
gamma.MG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between MG group and Control
gamma.RG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between RG group and Control
gamma.BG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between BG group and Control

# conditional distribution of the data (uses priors above)
for (i in 1:N){

# creating conditional mean to put into model statement 
score.hat[i] <- beta.0 + beta.exp4*exp4[i] + beta.enth*enthusiasm[i] + beta.MG*CvMG[i] + beta.RG*CvRG[i] + beta.BG*CvBG[i] + beta.exp4.MG*exp4[i]*CvMG[i] + beta.exp4.RG*exp4[i]*CvRG[i] + beta.exp4.BG*exp4[i]*CvBG[i] + beta.enth.MG*enthusiasm[i]*CvMG[i] + beta.enth.RG*enthusiasm[i]*CvRG[i] + beta.enth.BG*enthusiasm[i]*CvBG[i] + beta.exp4.enth*exp4[i]*enthusiasm[i] + beta.exp4.enth.MG*exp4[i]*enthusiasm[i]*CvMG[i] + beta.exp4.enth.RG*exp4[i]*enthusiasm[i]*CvRG[i] + beta.exp4.enth.BG*exp4[i]*enthusiasm[i]*CvBG[i]

# error values (for R^2)
error[i] = score[i]-score.hat[i]

varhat[i] = gamma.0 + gamma.MG*CvMG[i] + gamma.RG*CvRG[i] + gamma.BG*CvBG[i]
tau.e[i] <- exp(varhat[i])
# likelihood from model:
score[i] ~ dnorm(score.hat[i], tau.e[i])

}

dif.MG.RG <- beta.RG - beta.MG
dif.MG.BG <- beta.BG - beta.MG
dif.RG.BG <- beta.BG - beta.RG

difRPrec.MG.RG <- gamma.RG - gamma.MG
difRPrec.MG.BG <- gamma.BG - gamma.MG
difRPrec.RG.BG <- gamma.BG - gamma.RG

}
"

BayesianFullModel.Parameters = c("beta.0", 
                                 "beta.exp4",
                                 "beta.enth",
                                 "beta.MG",
                                 "beta.RG",
                                 "beta.BG",
                                 "beta.exp4.MG",
                                 "beta.exp4.RG",
                                 "beta.exp4.BG",
                                 "beta.enth.MG",
                                 "beta.enth.RG",
                                 "beta.enth.BG",
                                 "beta.exp4.enth",
                                 "beta.exp4.enth.MG",
                                 "beta.exp4.enth.RG",
                                 "beta.exp4.enth.BG",
                                 "dif.MG.RG",
                                 "dif.MG.BG",
                                 "dif.RG.BG",
                                 "difRPrec.MG.RG",
                                 "difRPrec.MG.BG",
                                 "difRPrec.RG.BG",
                                 "gamma.0", 
                                 "gamma.MG", 
                                 "gamma.RG", 
                                 "gamma.BG")


# in order to have reproducable Markov chains, we need to initialize the random number generator and seed values
BayesianFullModel.Seed = 12345678

#   Note: here, the random number seed cannot be the same per seed or the chains will be the same
RNGname = c("Wichmann-Hill","Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")
if (RNGname[1] %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                      "Super-Duper", "Mersenne-Twister")) {
  RNGname[1] <- paste("base::", RNGname[1], sep = "")
}

BayesianFullModel.Init.Values <- vector("list", n.chains)
for (i in 1:n.chains) {
  BayesianFullModel.Init.Values[[i]]$.RNG.name <- RNGname[1]
  BayesianFullModel.Init.Values[[i]]$.RNG.seed <- BayesianFullModel.Seed + i
}

load.module("glm")
load.module("dic")


# Submit model to JAGS for compiling. We'll turn adaptation off to control it in the next line of syntax:
BayesianFullModel.JAGS = jags.model(file = textConnection(BayesianFullModel.Syntax), data = BayesianFullModel.JAGS.Data, n.adapt = 0, n.chains = 4, inits = BayesianFullModel.Init.Values)

# adapting the model: No adaptation needed as Gibbs Sampling is being used -- but good to try in case algorithm isn't Gibbs
adapt(object = BayesianFullModel.JAGS, n.iter = 1000)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
BayesianFullModel.Samples = coda.samples(model = BayesianFullModel.JAGS, variable.names = BayesianFullModel.Parameters, n.iter = 10000, thin = 1)

# Remove burnin from samples:
BayesianFullModel.Posterior = window(x = BayesianFullModel.Samples, start = 5001, end = 10000)

hdi(BayesianFullModel.Posterior)

#plot(BayesianFullModel.Posterior)

#Calculates DIC value for this model
dic1<-dic.samples(BayesianFullModel.JAGS, n.iter = 1000)

summary(BayesianFullModel.Posterior)
hdi(BayesianFullModel.Parameters$beta.exp4.enth.MG)


#####Begin fit of Second Model#####

FullModel2 = lm(score~exp4+enthusiasm+CvMG+CvRG+CvBG+exp4:CvMG+exp4:CvRG+exp4:CvBG+enthusiasm:CvMG+enthusiasm:CvRG+enthusiasm:CvBG, data = BY02Data)

# take data frame and create model matrix -- which gives us the columns of X:
BayesianFullModel2.Matrix = model.matrix(object = FullModel2)

# append with weight and then select only columns needed for model:
BayesianFullModel2.Matrix = data.frame(cbind(BY02Data$score, BayesianFullModel2.Matrix))
names(BayesianFullModel2.Matrix) = c("score", "intercept", "exp4", "enthusiasm", "CvMG", "CvRG", "CvBG", "exp4*CvMG", "exp4*CvRG", "exp4*CvBG", "enthusiasm*CvMG", "enthusiasm*CvRG", "enthusiasm*CvBG")
BayesianFullModel2.Matrix = BayesianFullModel2.Matrix[c("score", "exp4", "enthusiasm", "CvMG", "CvRG", "CvBG")]

# add data and N to list object that we will pass to JAGS
BayesianFullModel2.JAGS.Data = list(N = nrow(BY02Data),
                                   score = BayesianFullModel.Matrix$score,
                                   exp4 = BayesianFullModel.Matrix$exp4,
                                   enthusiasm = BayesianFullModel.Matrix$enthusiasm,
                                   CvMG = BayesianFullModel.Matrix$CvMG,
                                   CvRG = BayesianFullModel.Matrix$CvRG,
                                   CvBG = BayesianFullModel.Matrix$CvBG,
                                   betaMean.0 = betaMean.0,
                                   betaTau.0 = betaTau.0)


BayesianFullModel2.Syntax = "

model{


# prior distributions: 
#   Note that terms on the left are model parameter names we create here
beta.0             ~ dnorm(betaMean.0, betaTau.0) # prior for the intercept
beta.exp4          ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional slope of exp4
beta.enth          ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional slope of enthusiasm
beta.MG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between MG group and Control
beta.RG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between RG group and Control
beta.BG            ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between BG group and Control
beta.exp4.MG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and MG group
beta.exp4.RG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and RG group
beta.exp4.BG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of exp4 and BG group
beta.enth.MG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and MG group
beta.enth.RG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and RG group
beta.enth.BG       ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of enthusiasm and BG group
gamma.0            ~ dnorm(betaMean.0, betaTau.0) #prior for the interaction of residual precision
gamma.MG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between MG group and Control
gamma.RG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between RG group and Control
gamma.BG           ~ dnorm(betaMean.0, betaTau.0) #prior for the conditional mean difference between BG group and Control

# conditional distribution of the data (uses priors above)
for (i in 1:N){

# creating conditional mean to put into model statement 
score.hat[i] <- beta.0 + beta.exp4*exp4[i] + beta.enth*enthusiasm[i] + beta.MG*CvMG[i] + beta.RG*CvRG[i] + beta.BG*CvBG[i] + beta.exp4.MG*exp4[i]*CvMG[i] + beta.exp4.RG*exp4[i]*CvRG[i] + beta.exp4.BG*exp4[i]*CvBG[i] + beta.enth.MG*enthusiasm[i]*CvMG[i] + beta.enth.RG*enthusiasm[i]*CvRG[i] + beta.enth.BG*enthusiasm[i]*CvBG[i]

# error values (for R^2)
error[i] = score[i]-score.hat[i]

varhat[i] = gamma.0 + gamma.MG*CvMG[i] + gamma.RG*CvRG[i] + gamma.BG*CvBG[i]
tau.e[i] <- exp(varhat[i])
# likelihood from model:
score[i] ~ dnorm(score.hat[i], tau.e[i])

}

dif.MG.RG <- beta.RG - beta.MG
dif.MG.BG <- beta.BG - beta.MG
dif.RG.BG <- beta.BG - beta.RG

difRPrec.MG.RG <- gamma.RG - gamma.MG
difRPrec.MG.BG <- gamma.BG - gamma.MG
difRPrec.RG.BG <- gamma.BG - gamma.RG

}
"

BayesianFullModel2.Parameters = c("beta.0", 
                                 "beta.exp4",
                                 "beta.enth",
                                 "beta.MG",
                                 "beta.RG",
                                 "beta.BG",
                                 "beta.exp4.MG",
                                 "beta.exp4.RG",
                                 "beta.exp4.BG",
                                 "beta.enth.MG",
                                 "beta.enth.RG",
                                 "beta.enth.BG",
                                 "dif.MG.RG",
                                 "dif.MG.BG",
                                 "dif.RG.BG",
                                 "difRPrec.MG.RG",
                                 "difRPrec.MG.BG",
                                 "difRPrec.RG.BG",
                                 "gamma.0", 
                                 "gamma.MG", 
                                 "gamma.RG", 
                                 "gamma.BG")


# in order to have reproducable Markov chains, we need to initialize the random number generator and seed values
BayesianFullModel2.Seed = 12345678

#   Note: here, the random number seed cannot be the same per seed or the chains will be the same
RNGname = c("Wichmann-Hill","Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")
if (RNGname[1] %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                      "Super-Duper", "Mersenne-Twister")) {
  RNGname[1] <- paste("base::", RNGname[1], sep = "")
}

BayesianFullModel2.Init.Values <- vector("list", n.chains)
for (i in 1:n.chains) {
  BayesianFullModel2.Init.Values[[i]]$.RNG.name <- RNGname[1]
  BayesianFullModel2.Init.Values[[i]]$.RNG.seed <- BayesianFullModel2.Seed + i
}

load.module("glm")
load.module("dic")


# Submit model to JAGS for compiling. We'll turn adaptation off to control it in the next line of syntax:
BayesianFullModel2.JAGS = jags.model(file = textConnection(BayesianFullModel2.Syntax), data = BayesianFullModel2.JAGS.Data, n.adapt = 0, n.chains = 4, inits = BayesianFullModel2.Init.Values)

# adapting the model: No adaptation needed as Gibbs Sampling is being used -- but good to try in case algorithm isn't Gibbs
adapt(object = BayesianFullModel2.JAGS, n.iter = 1000)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
BayesianFullModel2.Samples = coda.samples(model = BayesianFullModel2.JAGS, variable.names = BayesianFullModel2.Parameters, n.iter = 10000, thin = 1)

# Remove burnin from samples:
BayesianFullModel2.Posterior = window(x = BayesianFullModel2.Samples, start = 5001, end = 10000)

#plots the posterior
pdf("Trace Plots.pdf")
plot(BayesianFullModel2.Posterior)
dev.off()

summary(BayesianFullModel2.Posterior)

hdi(BayesianFullModel2.Posterior)

#calculates DIC
dic2<-dic.samples(BayesianFullModel2.JAGS, n.iter = 1000)

#calculates differenc of DIC values from model 1 and model 2
diffdic(dic1,dic2)

#Gelman Rubin Diagnostic
gelman.diag(BayesianFullModel2.Posterior, multivariate = FALSE)

#produces autocorrelation information
autocorr.diag(BayesianFullModel2.Posterior)

#produces 0.95 HDI for parameters
hdi(BayesianFullModel2.Posterior)
