
##############################################
### Load data and produce the x variable
##############################################

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


res <- LoadAndGenerate("AntMAN/tests/sthree.jpg")
RunKmeans(res$x,res$dim)


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

##############################################
### mat stupid
##############################################



matstupid_n=1000
matstupid=matrix(0,matstupid_n,3)
matstupid[,1]<- sample(c(0.1,.8),matstupid_n,replace=TRUE)
matstupid[,2]<- sample(c(0.3,.7),matstupid_n,replace=TRUE)
matstupid[,3]<- sample(c(0.2,.5),matstupid_n,replace=TRUE)
scatter3D(x=matstupid[,1],y=matstupid[,2],z=matstupid[,3])

##############################################
### random three dim
##############################################


MU <- matrix(nrow=3,ncol=3)

MU[1,] <- c(0,0,0)
MU[2,] <- c(-3,-3,-3)
MU[3,] <- c(4,4,4)


sig1 <- c(1,1,1)
rho1 <- 0
Sig1 <- diag(1,3)


sig2 <- c(1,1)
rho2 <- -0.7
Sig2 <- matrix(c(sig2[1]^2,rho2*sig2[1]*sig2[2], rho2*sig2[1]*sig2[2],sig2[2]^2),byrow=T,nrow=2) 
Sig2<- matrix(c(1,-0.7,0,-0.7,1,-0.7,0,-0.7,1),3,3)
Sig2

sig3 <- c(1,1,1)
rho3 <- -0.3
Sig3 <- matrix(c(sig3[1]^2,rho3*sig3[1]*sig3[2], rho3*sig3[1]*sig3[2],sig3[2]^2),byrow=T,nrow=2) 

Sig3<- matrix(c(1,-0.3,-0.3,-0.3,1,-0.3,-0.3,-0.3,1),3,3)
Sig3

SIG <- array(0,dim=c(3,3,3))
SIG[1,,] <- Sig1
SIG[2,,] <- Sig2
SIG[3,,] <- Sig3


library(mvtnorm) ##  TODO Why do I need this !!!!
demo_multivariate_normal <-AM_sample_multinorm(n = 1000 ,d = 3,c(0.3,0.3,0.4),MU,SIG)
y_mvn  <- demo_multivariate_normal$y
ci_mvn <- demo_multivariate_normal$ci

hist(y_mvn,freq=F,nclass=15,col=colors()[4])

scatter3D(x=y_mvn[,1],y=y_mvn[,2],z=y_mvn[,3])

##############################################################################
### PREPARE THE GIBBS for multivariate Normal mixture with poisson gamma priors
##############################################################################
