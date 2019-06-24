
##############################################
### Load packages
##############################################
.libPaths(c("./AntMAN.Rinstall", .libPaths()))

library("AntMAN")


rm(list=ls())
set.seed(123)

##############################################
### Functions
##############################################

  


DrawResult <- function(mat,imgDm,clusters) {
  clusters = clusters / max(clusters)
  img_seg <- array(dim=imgDm)
  fit$CI[-1]
  img_seg[,,1] <- matrix(clusters,nrow=imgDm[1],ncol=imgDm[2])
  img_seg[,,2] <- matrix(clusters,nrow=imgDm[1],ncol=imgDm[2])
  img_seg[,,3] <- matrix(clusters,nrow=imgDm[1],ncol=imgDm[2])
  plot(1:2, type='n')
  rasterImage(img_seg, 1.01, 1.01, 1.99, 1.99)
  
  
  
  
}

data(brain)
imgDm <- brain$dim
x = brain$pic
bn <- imgDm[1] * imgDm[2]
mat <- matrix(0,bn,3)
mat[,1 ] <-x$R
mat[,2 ] <-x$G
mat[,3 ] <-x$B
### scatter3D(x=mat[,1],y=mat[,2],z=mat[,3])


mixture_mvn_params = AM_multinorm_mix_hyperparams   (mu0=c(0,0,0),ka0=.1,nu0=5,Lam0=0.1*diag(3))

mcmc_params        = AM_mcmc_parameters(niter=1000, burnin=0, thin=5, verbose=3, output = c("CI","K","M"))
components_prior   = AM_mix_components_prior_pois (init=5,a=10,b=2) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
  init_K = 1,
  y = mat, 
  mix_kernel_hyperparams = mixture_mvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

fres = fit$CI[[length(fit$CI)]]
DrawResult(mat,imgDm,fres)

summary (fit)
plot (fit)

