# preliminary to extract the data from the poLCA package, to be cancelled
# install.packages("poLCA")
# library(poLCA)
# data(carcinoma)
# y <- carcinoma-1
# head(y) 
# write.matrix(file="carcinoma.txt",y)


### Read the data. This dataset should be part of AntMAN
carc <- as.matrix(read.table("carcinoma.txt",header=T))

# Basic latent class modeling with the carcinoma data
# The carcinoma data from Agresti (2002, 542) consist
# of seven dichotomous variables that represent the 
# ratings by seven pathologists of 118 slides on the presence or 
# absence of carcinoma.

# The purpose of studying these data is to model "interobserver agreement"
# by examining how subjects might be divided into groups depending upon 
# the consistency of their diagnoses.

# It is straightforward to replicate Agresti's published results 
# (Agresti 2002, 543) using the
# series of commands:






### Da Agresti 
#Seven pathologists classified each of 118 slides on the presence or absence of
#carcinoma in the uterine cervix. For modeling interobserver agreement, the
#conditional independence assumption of the latent class model is often
#plausible. With a blind rating scheme, ratings of a given subject or unit by
#different pathologists are independent. If subjects having true rating in a
#given category are relatively homogeneous, then ratings by different patholo-
#gists may be nearly independent within a given true rating class. Thus, one
#might posit a latent class model with q s 2 classes, one for subjects whose
#true rating is positive and one for subjects whose true rating is negative. This



# model expresses the 2 7 joint distribution of the seven ratings 
# as a mixture of
# two 2 7 distributions, one for each true rating class.
# Table 13.2 shows results of fitting some latent class models Žincluding a
# mixture model studied in Section 13.2.4.. Because the observed table is
# sparse, the deviance is mainly useful for comparing models. This is an
# informal comparison, though, since the chi-squared distribution does not
# apply for comparing deviances of models with different numbers of latent
# classes. A model with q classes is a special case of a model with q* ) q
# classes in which P Ž Z s z . s 0 for z ) q and hence falls on the boundary of
# the parameter space. Ordinary chi-squared likelihood-ratio tests require
# parameters to fall in the interior of the parameter space 
# Ži.e., 0 - P Ž Z s z .
# - 1 for z s 1, . . . , q*..










y <- carc

## Just verify that we have loaded the right data
is(y)
head(y)

### Data quantities
n <- dim(y)[1] ## The number of observation
d <- dim(y)[2] ## The dimension of the data



## Load AntMan
#library("AntMan", lib.loc = "~/Rcpp/finite_mixture/AntMan.install")
cartella <- getwd()
setwd("/home/bruffo/Desktop/github/")
library("AntMAN", lib.loc = "AntMAN.Rinstall")
setwd(cartella)


######
## Set the Gibbs Parameter
mcmc_params        = AM_mcmc_parameters(niter=5000, burnin=1000, thin=10, verbose=1, output = c("CI","K","M","TAU"))

## We are going to use a Bayesian Latent Class analysis, i.e. a mixture of Multivariate Bernoulli
## for the first analysis I will define independent beta prior with (1,1) parameters (i.e. uniform)
## for each component
mixture_mvb_params <- AM_multiber_mix_hyperparams(a0=rep(1,d),b0= rep(1,d))
components_prior   <-  AM_mix_components_prior_pois (init=5,a=10,b=2) 
#components_prior   <-  AM_mix_components_prior_dirac (Mstar=30)
weights_prior      <-  AM_mix_weights_prior_gamma(init=2, a=1, b=1)



dim(y)[1]
### Let's start with all the data in separate cluster (i.e. n groups each of size 1)
init_ci_mvb <-  0:(n-1) ## kmeans(y,centers=4)$cluster ### 0:(n-1)
fit <- AM_mcmc_fit(
  init_K = 1,
  y = as.matrix(y), 
  mix_kernel_hyperparams = mixture_mvb_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary.AM_mcmc_fitness_result(fit)

plot(fit$K)
                             


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

