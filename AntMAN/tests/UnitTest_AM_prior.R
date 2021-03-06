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
### Run
##############################################

prior_K_de = AM_prior_K_Delta(1000,0.1,30)

##############################################
### Verify
##############################################

stopifnot(class(prior_K_de) == "AM_prior")
