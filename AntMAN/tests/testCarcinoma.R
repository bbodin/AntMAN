pdf("Carcinoma.pdf")

.libPaths(c("./AntMAN.Rinstall", .libPaths()))

library("AntMAN")

data(carcinoma)

### Data quantities
y <- carcinoma
n <- dim(y)[1] ## The number of observation
d <- dim(y)[2] ## The dimension of the data
mcmc_params        = AM_mcmc_parameters(niter=5000, burnin=1000, thin=10, verbose=1, output = c("CI","K","M","TAU"))

## We are going to use a Bayesian Latent Class analysis, i.e. a mixture of Multivariate Bernoulli
## for the first analysis I will define independent beta prior with (1,1) parameters (i.e. uniform)
## for each component
mixture_mvb_params <- AM_multiber_mix_hyperparams(a0=rep(1,d),b0= rep(1,d))
components_prior   <-  AM_mix_components_prior_pois (init=5,a=10,b=2) 
#components_prior   <-  AM_mix_components_prior_dirac (Mstar=30)
weights_prior      <-  AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
  init_K = 1,
  y = as.matrix(y), 
  mix_kernel_hyperparams = mixture_mvb_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary(fit)
plot(fit)
readline(prompt="Press [enter] to continue")
print(fit$TAU)



ci <- t(do.call(cbind,fit$CI))

dim(ci)

library('sdols')
clu <- dlso(ci)

clu

## estimated number of clusters
hatk <- length(unique(clu))
est <- matrix(ncol=dim(y)[2],nrow=hatk)

#A method for determining the size of the latent classes is to assign
# each observation to a latent class on an individual basis according to its model posterior
# class membership probability. Values using this technique are reported directly below the
#estimated mixing proportions. Congruence between these two sets of population shares often
#indicates a good fit of the model to the data.



par(mfrow=c(1,3))
for(j in 1:hatk){
  est[j,] <- apply(y[clu==j,],2,mean)
  plot(apply(y[clu==j,],2,mean),type="h",xaxt="n",xlim=c(0.8,6.2),ylim=c(0,1),col=j,lwd=2,ylab="")
  axis(1, at=1:7, labels=colnames(y))
}


est



table(clu)/length(clu)
#The figure shows  As Agresti describes, the three estimated latent classes clearly correspond to
#a pair of classes that are consistently rated negative (37%) or positive (44%), plus a third
#2problematic" class representing 18% of the population. 
colnames(est) <- names(y)

est[2,]
#In that class, pathologists B, E, and  G tend to diagnose positive; C, D, and F tend to diagnose
# negative; and A is about 50/50.

##DEVO FARLO
# First, the estimated classcconditional response probabilities p_rjk are reported for pathologists 
# A through G, with each row corresponding to a latent class, and each column corresponding to a diagnosis; negative in the first column, and positive in the second. 

#Thus, for example, a slide belonging to the first (“negative”) class has a 94% chance of being rated free from carcinoma by rater A, an 86% chance of the same from rater B, an 100% chance from rater C, and so forth.


#Next, the output provides the estimated mixing proportions p̂r corresponding to the share
#of observations belonging to each latent class. These are the same values that appear in
# Figure 1. 



names(fit)

fit$TAU
