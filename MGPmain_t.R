if (!require("Rcpp", character.only = TRUE)) {
  install.packages("Rcpp", dependencies = TRUE)
}

# set working directory as the current location of this MGPmain_t.R file.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(Rcpp)
sourceCpp("mgp_funcs_t.cpp")
sourceCpp("mgp_funcs.cpp")
source("MGPprior_t.R")
source("MGPprior.R")

# Simulation settings ########################################################################
num.vars = 500 # number of variables (p)
num.factors = 15 # number of latent factors (k)
n = 200 # number of observations (n)
nu = 3 # True degrees of freedom of t-distribution

# Generate simulation data ###################################################################
datagen = generateData(num.vars = num.vars, num.factors = num.factors, n = n, 
                       student_t = TRUE, nu = nu, sparse = TRUE)
X = datagen[[1]]
true_cov = datagen[[2]]

sum(true_cov == 0)/(num.vars*num.vars) # ratio of zero entries in true covariance.

#################################################################################################
# t
nuset = 3 # prespecified degrees of freedom of t-likelihood of factor model.
start = Sys.time()
result = linearMGSP_t(as.matrix(X), nrun = 20000, burn = 5000, thin = 1, nu = nuset, prop = 0.7, 
                      epsilon = 0.01, kinit = 50, output = c("covMean","numFactors","factSamples"),
                      epseta = 0.015)
Sys.time() - start

#################################################################################################
# Normal
start = Sys.time()
result2 = linearMGSP(as.matrix(X), nrun = 20000, burn = 5000, thin = 1, prop = 0.7, 
                     epsilon = 0.01, kinit = 50, output = c("covMean","numFactors"))
Sys.time() - start

#################################################################################################
mean((true_cov - result$covMean)^2)  # Mean Squared Error of t factor model
mean(abs(true_cov - result$covMean))  # Mean Absolute Deviation of t factor model
max(abs(true_cov - result$covMean))  # Maximum Absolute Deviation of t factor model
print("################################################")
mean((true_cov - result2$covMean)^2)  # Mean Squared Error of normal factor model
mean(abs(true_cov - result2$covMean))  # Mean Absolute Deviation of normal factor model
max(abs(true_cov - result2$covMean))  # Maximum Absolute Deviation of normal factor model
########################################################
########################################################
norm(result$covMean - true_cov,"1") # one norm (maximum absolute column sum) of t factor model
norm(result$covMean - true_cov,"2") # two norm (maximum sigular value) of t factor model
norm(result$covMean - true_cov,"i") # infinity norm (maximum absolute row sum) of t factor model
norm(result$covMean - true_cov,"F") # Frobenius norm (L2norm of vec(cov-cov)) of t factor model
norm(result$covMean - true_cov,"M") # Maximum modulus (Maximum Absolute Deviation) of t factor model
########################################################
norm(result2$covMean - true_cov,"1") # one norm (maximum absolute column sum) of normal factor model
norm(result2$covMean - true_cov,"2") # two norm (maximum sigular value) of normal factor model
norm(result2$covMean - true_cov,"i") # infinity norm (maximum absolute row sum) of normal factor model
norm(result2$covMean - true_cov,"F") # Frobenius norm (L2norm of vec(cov-cov)) of normal factor model
norm(result2$covMean - true_cov,"M") # Maximum modulus (Maximum Absolute Deviation) of normal factor model
########################################################
# 2.5 and 97.5 percentile of estimates for true zero entries.
quantile(as.vector(result$covMean[true_cov == 0]), c(0.025,0.975))
quantile(as.vector(result2$covMean[true_cov == 0]), c(0.025,0.975))

# mean number of latent factors
mean(result$numFacts) # t factor model
mean(result2$numFacts) # normal factor model

# few quantities for MCMC diagnostic of robust sparse Bayesian infinite factor model.
par(mfrow = c(5,1), mar=c(2,4,1,3))
plot(result$numFacts,type='l',main = "MGP - t", ylab='# of factors')
plot(apply(result$gammaSamps,2,mean),type='l', ylab='gam')
plot(log2(result$treecount),type='l', ylab = 'log_2(treecount)')
plot(result$a1Samps,type='l', ylab='a1')
plot(result$a2Samps,type='l', ylab='a2')

# Plot of entries of true covariance with entries of estimated covariance.
par(mfrow = c(1,2))
lim = max(max(abs(as.vector(result$covMean))),max(abs(as.vector(result2$covMean))))
xylim = c(-lim, lim)
plot(as.vector(true_cov), as.vector(result$covMean), xlab="true value", ylab="estimated value", 
     xlim = xylim, ylim = xylim, main = "MGP - student t",pch = 16, col='#0000000A')
abline(0,1)
abline(v=0)
abline(h=0)

plot(as.vector(true_cov), as.vector(result2$covMean), xlab="true value", ylab="estimated value", 
     xlim = xylim, ylim = xylim, main = "MGP - normal",pch = 16, col='#0000000A')
abline(0,1)
abline(v=0)
abline(h=0)
