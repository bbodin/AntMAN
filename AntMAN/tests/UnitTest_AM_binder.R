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

eam = AM_clustering(res$fit)
cluster = AM_salso(eam, "binder")

# binder_result = AM_binder(res$fit)
# cluster = binder_result[["Labels"]]
# coclustering_probability = AM_coclustering(res$fit)
#binder_result = AM_binder(res$fit , with_coclustering_probability=TRUE)
# cluster                    = binder_result[["clustering"]]
# coclustering_probability   = binder_result[["coclustering_probability"]]

##############################################
### Verify
##############################################

stopifnot(is.vector(cluster))
stopifnot(is.matrix(eam))