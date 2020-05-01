#######################################################################################
###############
############### AntMAN Package : Tests and Examples
###############
###############
#######################################################################################



library("AntMAN")

prior_K_de = AM_prior_K_Delta(1000,0.1,30)
#assert( "Ensure the return type of AM_prior_K_Delta is valid.", class(prior_K_de) == "AM_prior_K")