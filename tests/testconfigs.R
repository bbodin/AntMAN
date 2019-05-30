library(testthat)
library(AntMan)

context("Configuration generation") 

test_that("AM_mcmc_parameters return a list a parameters", {
    expect_equal(AM_mcmc_parameters()$type,"AM_MCMC_PARAMETERS") 
}) 



