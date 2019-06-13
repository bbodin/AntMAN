
##############################################
### Load packages
##############################################
.libPaths(c("./AntMAN.Rinstall", .libPaths()))
library("jpeg")
library(plot3D)
library("AntMAN")


rm(list=ls())
set.seed(123)

##############################################
### Functions
##############################################

### Load image and generate dataset
LoadAndGenerate <- function(img_path) {
  
  # Import the image in R
  img <- readJPEG(img_path) # Read the image
  
  ## We can plot the image with the command rasterImage
  plot(1:2, type='n')
  rasterImage(img, 1.01, 1.01, 1.99, 1.99)
  
  # Obtain the dimension
  imgDm <- dim(img)
  
  # Assign RGB channels to data frame
  imgRGB <- data.frame(
    x = rep(1:imgDm[2], each = imgDm[1]),
    y = rep(imgDm[1]:1, imgDm[2]),
    R = as.vector(img[,,1]),
    G = as.vector(img[,,2]),
    B = as.vector(img[,,3])
  )
  
  x <- imgRGB[,c(3,4,5)]
  
  return (list(x = x, dim = imgDm))
}

RunKmeans <- function(x , imgDm , kClusters = 5) {

  
  ### cluster the pixels according the RGB channels
  kMeans <- kmeans(x, centers = kClusters)
  
  
  ### This are the k centroids that we have (the k colours with which we aim at representing the image)
  
  rgb(kMeans$centers)
  
  ### The image reconstructed with just k colors (the kcentroids)
  segmented <- kMeans$centers[kMeans$cluster,]
  
  ### we have to transform the data with three colums in an an array 
  ### of three matrices 
  img_seg <- array(dim=imgDm)
  
  img_seg[,,1] <- matrix(segmented[,1],nrow=imgDm[1],ncol=imgDm[2])
  img_seg[,,2] <- matrix(segmented[,2],nrow=imgDm[1],ncol=imgDm[2])
  img_seg[,,3] <- matrix(segmented[,3],nrow=imgDm[1],ncol=imgDm[2])
  plot(1:2, type='n')
  rasterImage(img_seg, 1.01, 1.01, 1.99, 1.99)
  
}


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

res <- LoadAndGenerate("AntMAN/tests/peppers100.jpg")
x <- res$x
imgDm <- res$dim
bn <- imgDm[1] * imgDm[2]
mat <- matrix(0,bn,3)
mat[,1 ] <-x$R
mat[,2 ] <-x$G
mat[,3 ] <-x$B

scatter3D(x=mat[,1],y=mat[,2],z=mat[,3])


mixture_mvn_params = AM_multinorm_mix_hyperparams   (mu0=c(0,0,0),ka0=.1,nu0=5,Lam0=0.1*diag(3))

mcmc_params        = AM_mcmc_parameters(niter=500, burnin=0, thin=5, verbose=1, output = c("CI","K","M"))
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

