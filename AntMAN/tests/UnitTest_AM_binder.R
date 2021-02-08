#######################################################################################
###############
############### AntMAN Package : Tests and Examples
###############
###############
#######################################################################################

##############################################
### Load the AntMan package
##############################################

library("AntMAN")
set.seed(123)

##############################################
### Prepare the data
##############################################

res = AM_demo_uvn_poi () 

##############################################
### Run
##############################################

binder_result = AM_salso(res$fit, with_coclustering_probability=TRUE)
#binder_result = AM_binder(res$fit , with_coclustering_probability=TRUE)

cluster                    = binder_result[["clustering"]]
coclustering_probability   = binder_result[["coclustering_probability"]]

##############################################
### Verify
##############################################

stopifnot(is.vector(cluster))
stopifnot(is.matrix(coclustering_probability))