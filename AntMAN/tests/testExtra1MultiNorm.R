###########################################################
# Same, but with hihg-dimensional data #
#*********************************************************#

rm(list=ls())

library("mvtnorm")
library("AntMAN")
getNamespaceVersion("AntMAN")

set.seed(123)

p <- 10
MU <- matrix(nrow=3,ncol=p)

MU[1,] <- rmvnorm(n=1, rep(0,p), diag(1,p))
MU[2,] <- MU[1,]
MU[3,] <- MU[1,]

Sig1 <- diag(0.5,p)
Sig2 <- diag(0.55,p)
Sig3 <- diag(0.45,p)

SIG <- array(0,dim=c(3,p,p))
SIG[1,,] <- Sig1
SIG[2,,] <- Sig2
SIG[3,,] <- Sig3


N <- 1000
demo_multivariate_normal <- AM_sample_multinorm(n = N,d = p,c(0.3,0.3,0.4),MU,SIG)
y_mvn  <- demo_multivariate_normal$y
ci_mvn <- demo_multivariate_normal$ci

hist(y_mvn,freq=FALSE,nclass=15,col=colors()[4])
plot(y_mvn,col=ci_mvn+1)



###
# Run Mixture model with Split/Merge algorithm from our function
y <- y_mvn


a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1
Lambda <- a2/b2
gamma_S <- a1/b1


#MCMC length
n_burn1 <- 0
n_burn2 <- 0
n_save <- 1000
thin <- 1
niter <- n_burn1 + n_burn2 + n_save * thin

mcmc_params        = AM_mcmc_parameters(niter=niter, burnin=n_burn1 + n_burn2, thin=thin, verbose=1, output = c("CI","K","M","H","Q"))
components_prior   = AM_mix_components_prior_pois(Lambda = Lambda)
weights_prior      = AM_mix_weights_prior_gamma(gamma = gamma_S)

#Initialize allocations
c_init <- kmeans(x = y, centers = p)
c_init <- c_init$cluster
# c_init <- c(1:N)
# c_init <- rep(1,N)

#Hyperparameters
mu0 <- rmvnorm(n=1, rep(-1,p), diag(1,p))
k0 <- 1
nu0 <- 2*(p + 2)
Lam0 <- diag(p)
# Lam0 <- (nu0 - p - 1) * Sig

mixture_mvn_params = AM_multinorm_mix_hyperparams(mu0=mu0,ka0=k0,nu0=nu0,Lam0=Lam0)


#Fit AntMAN (different hyperparameters!!!)
set.seed(321)
MCMC_output_AntMAN <- AM_mcmc_fit(
  initial_clustering = c_init,
  y = y,
  mix_kernel_hyperparams = mixture_mvn_params,
  mix_components_prior = components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary(MCMC_output_AntMAN)

