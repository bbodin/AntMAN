quit() ## Skip this test - too long
###########################################################
# Same, but with hihg-dimensional data #
#*********************************************************#

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

mcmc_params        = AM_mcmc_parameters(niter=niter, burnin=n_burn1 + n_burn2, thin=thin, verbose=0, output = c("ALL"))
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

mixture_mvn_params = AM_mix_hyperparams_multinorm(mu0=mu0,ka0=k0,nu0=nu0,Lam0=Lam0)


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


loglike_MCMC_iter = rep(0,niter)

for (it in 1:niter) {
  
  loglike_MCMC = 0
  c = MCMC_output_AntMAN$CI[[it]]
  S_m = MCMC_output_AntMAN$S[[it]]
  mu_star = MCMC_output_AntMAN$TAU[[it]]$mu
  Sigma_star = MCMC_output_AntMAN$TAU[[it]]$Sig
  M = MCMC_output_AntMAN$M[[it]]

  
  for (i in 1:N) {
    loglike_MCMC = loglike_MCMC + log(S_m[c[i]+1]) + dmvnorm(y[i,], mu_star[c[i]+1,], Sigma_star[,,c[i]+1], TRUE);
  }
  for (j in 1:M) {
    loglike_MCMC = loglike_MCMC + dmvnorm(mu_star[j,], mu0, Sigma_star[,,j]/k0, TRUE);
    loglike_MCMC = loglike_MCMC + log(diwish( Sigma_star[,,j], nu0, Lam0));
    loglike_MCMC = loglike_MCMC + dgamma(S_m[j], gamma_S, scale=1,log= TRUE);
  }
  
  loglike_MCMC = loglike_MCMC + dgamma(Lambda, a2, scale=1/b2, log=TRUE);
  loglike_MCMC = loglike_MCMC + dgamma(gamma_S, a1, scale=1/b1,log= TRUE);
  #loglike_MCMC = loglike_MCMC + dgamma(u, N, , rate=sum(S_m));
  
  loglike_MCMC_iter[it] = loglike_MCMC
  
  print(loglike_MCMC)
}
