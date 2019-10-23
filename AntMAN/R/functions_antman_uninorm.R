### This function takes as input the desired values of
### the marginal mean of mu = Emu
### the marginal variance of mu = Vmu
### the marginal mean of sig2 = Esig2 
### the marginal variance of sig2 = Vsig2
###
### The values returns a list with 
### hyperparameters of a Normal-Inverse-Gamma
### with the desired value of the marginal moments.
from_statitstics2hyperparam <- function(Emu,Vmu,Esig2,Vsig2){

	m0    = Emu
	nu0   = 2*(Esig2)^2/Vsig2+4
	sig02 = Esig2*(nu0-2)/nu0
	k0    = sig02/Vmu * nu0/(nu0-2)

	return(list(m0=m0,nu0=nu0,sig02=sig02,k0=k0))
}


### this function compute the hyperparameters of an Normal-Inverse-Gamma 
### distribution using an empirical Bayes approach.
### The data y
### a positive value scEmu (default =1) such that marginally E(mu)=(sample variance)*scEmu
### a positive value scEsig2 (default=3) such that marginally E(sig2)=(sample variance)*scEsig2
### The coefficient of variation of sig2 (default =3), CVsig2
### 

emp_bayes_uninorm <- function(y,scEmu=1,scEsig2=3,CVsig2=3){
	n <- length(y)   ### sample size
	bary <- mean(y)  ### sample mean
	s2y <- var(y)    ### sample variance

	Emu <- bary
	Vmu <- s2y*scEmu
	Esig2 <- s2y/scEsig2
	Vsig2 <- CVsig2^2*Esig2^2

	return(from_statitstics2hyperparam(Emu,Vmu,Esig2,Vsig2))

}







## This function  produce nice representation of the posterior chain 
## of univariate hyperparameter


univariate_plot <- function(chain,  main_trace="Trace plot", main_auto="Autocorrelation",main_hist="Histogram",col_hist="gray",col_dens="blue",col_quant="red",nclass="fd",lag_acf=30){
### A nice plot to  synthesize our analysis on the parameter theta!
### Please try to understand the code
#x11()

	if(length(unique(chain))==1){
		warning("The chain is constant")
		stop();
	}




## Representation of the posterior chain of  theta
#Divide the plot device in three sub-graph regions
#two square on the upper and a rectangle on the bottom
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
#trace-plot of the posterior chain
plot(chain,type="l",main=main_trace)
# autocorrelation plot
acf(chain,main=main_auto,lag.max=lag_acf)


#Histogram
hist(chain,nclass=nclass,freq=FALSE,main=main_hist,col=col_hist) 

## Overlap the kernel-density 
lines(density(chain),col=col_dens,lwd=2)

## Posterior credible interval of beta0
## quantile(chain,prob=c(0.05,0.5,0.95))


## We Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.05)),col=col_quant,lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col=col_quant,lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.95)),col=col_quant,lty=2,lwd=2)

## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c(col_quant,col_quant,col_dens),lty=c(1,2,1))


layout(matrix(c(1,1,1,1),1,1,byrow=TRUE))
}








###### Predictive:







