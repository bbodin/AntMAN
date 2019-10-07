
#################################################################################
##### Raffaele Functions
#################################################################################

#' Compute the logarithm of the absolute value of the  generalized Sriling number of second Kind (mi pare) See charambeloides, using a recursive formula Devo vedere la formula
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param gamma    A positive real number \code{gamma} 
#'
#' @return A vector of length n, reporting the values \code{C(gamma,n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' dd= AM_compute_stirling_ricor_abs(11,10)
#' print(dd)
AM_compute_stirling_ricor_abs <- function (n,gamma) {
	return(compute_stirling_ricor_abs(n,gamma));
}


#' Compute ...
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param gamma    A positive real number \code{gamma} 
#'
#' @return A vector of length n, reporting the values ... for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' dd= AM_compute_stirling_ricor_log(11,10)
#' print(dd)
AM_compute_stirling_ricor_log<- function (n,gamma) {
	return(compute_stirling_ricor_log(n,gamma));
}

#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Shifted Poisson of parameter \code{Lambda}. See Section 9.1.1 of the Paper Argiento de Iorio 2019.
#' 
#' There are no default values.
#'
#' @param n        The sample size
#' @param Lambda   The \code{Lambda} parameter of the Poisson
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=1000
#' Lam=100
#' gam=0.5
#' vnk= AM_VnkPoisson(n,Lam,gam)
#' stir= AM_compute_stirling_ricor_log(n, gam)
#' plot(exp(vnk+stir))
#' sum(exp(vnk+stir ))

AM_VnkPoisson <- function (n,Lambda,gamma) {
	return(VnkPoisson(n,Lambda,gamma));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Shifted Poisson of parameter \code{Lambda}. See Section 9.1.1 of Argiento de Iorio (2019) for more details.
#' 
#' There are no default values.
#'
#' @param n        The sample size
#' @param Lambda   The \code{Lambda} parameter of the Poisson
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values of the prior on the number of clusters induced by the prior on \code{M} and \code{w}, i.e. \code{p^*_k} for \code{k=1,...,n}. See Section 9.1.1 of Argiento de Iorio (2019) for more details.
#'
#' @keywords prior number of clusters
#'
#' @export
#' 
#' @examples
#' n <- 82
#' Lambda <- 10
#' gam_po <- 0.1550195
#' prior_K_po <-  AM_prior_K_Pois(n,gam_po,Lambda)
#' plot(1:n, prior_K_po, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)")

AM_prior_K_Pois <- function (n,gamma,Lambda) {
	return(prior_K_Pois(n,gamma,Lambda));
}


#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet is Negative Binomial with parameter \code{r} and \code{p}with  mean is mu =1+ r*p/(1-p) TODO: CHECK THIS FORMULA!!!. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param gam    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=1000
#' r=1000
#' p=0.5
#' gam=0.5
#' vnk= AM_VnkNegBin(n,r,p,gam);
#' stir= AM_compute_stirling_ricor_log(n, gam)
#' plot(exp(vnk+stir+(1:n)*log(gam)))
#' sum(exp(vnk+stir))

AM_VnkNegBin <- function (n,r,p,gam) {
	return(VnkNegBin(n,r,p,gam));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) ) when the prior on the number of componet  is Negative Binomial with parameter \code{r>0} and \code{0<p<1}, with  mean is mu =1+ r*p/(1-p) TODO: CHECK THIS FORMULA!!!. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details. 
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param gamma  The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 50
#' gamma <- 1
#' r <- 0.1
#' p <- 0.91
#' gam_nb <- 0.2381641
#' prior_K_nb <-  AM_prior_K_NegBin(n,gam_nb,r,p)
#' plot(1:n,prior_K_nb, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)")
#' lines(1:n,prior_K_nb,type="h",lwd=2)

AM_prior_K_NegBin <- function (n,gamma, r, p){
	return(prior_K_NegBin(n,gamma, r, p));
}


#' Compute the value V(n,k), needed to caclulate the eppf of a Finite Dirichlet process when the prior on the component-weigts
#'  of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights are distributed as Gamma(\eqn{\gamma},1) )
#'   when the number of component are fixed to \code{M^*}, i.e. a Dirac prior assigning mass only to \code{M^*} is assumed. 
#'   See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details.
#' 
#' There are no default values.
#'
#' @param n      The sample size
#' @param Mstar  The number of component of the mixture 
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n=200
#' Mstar=100
#' gam=0.5
#' vvv=AM_VnkDelta(n,Mstar,gam);
#' stir= AM_compute_stirling_ricor_log(n, gam)
#' stir
#' plot(exp(vvv+stir) )
#' sum(exp(vvv+stir ))

AM_VnkDelta <- function (n,Mstar,gamma) {
	return(VnkDelta(n,Mstar,gamma));
}

#' This function compute the prior on the number of cluster, i.e. occupied component of the mixutre for a Finite Dirichlet process 
#' when the prior on the component-weigts of the mixture is a Dirichlet with parameter \code{gamma} (i.e. when unnormailized weights 
#' are distributed as Gamma(\eqn{\gamma},1) ) when the number of component are fixed to \code{M^*}, i.e. a Dirac prior assigning mass
#'  only to \code{M^*} is assumed. See Section 9.1.1 of the Paper Argiento de Iorio 2019 for more details.#' There are no default values.
#'
#' @param n        The sample size
#' @param Mstar    The number of component of the mixture 
#' @param gamma    The \code{gamma} parameter of the Dirichlet 
#'
#' @return A vector of length n, reporting the values \code{V(n,k)} for \code{k=1,...,n}
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' gam_de <- 0.1743555
#' Mstar <- 12
#' prior_K_de <- AM_prior_K_Delta(n,gam_de, Mstar)
#' plot(1:n, prior_K_de, type = "n", bty = "l", xlab = "k", ylab = "P(K=k)")

AM_prior_K_Delta <- function (n,gamma,Mstar){
	return( prior_K_Delta(n,gamma,Mstar));
}


#' Once specified a fixed value of components \code{M^*} this function  adopt a  \emph{bisection method} to find the value of \code{gamma} 
#' such that the induced distribution on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use
#'  a bisection method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} and
#'  an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.
#'  A default value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31), the function 
#'  stops warning that convergence has not bee reached.
#'
#' @param n             sample size
#' @param Mstar         number of component of the mixture 
#' @param Kstar         mean number of cluster the user want to specify
#' @param gam_min       lower bound of the interval in which \code{gamma} should be lie (default 1e-4)
#' @param gam_max       upper bound of the interval in which \code{gamma} should lie (default 10)
#' @param tolerance     tolerance for the method
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' Mstar <- 12
#' gam_de <- AM_find_gamma_Delta(n,Mstar,Kstar=6, gam_min=1e-4,gam_max=10, tolerance=0.1)
#' prior_K_de <-  AM_prior_K_Delta(n,gam_de,Mstar)
#' prior_K_de\%*\%1:n

AM_find_gamma_Delta <- function (n,Mstar,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1) {
	return(find_gamma_Delta(n,Mstar,Kstar, gam_min,gam_max, tolerance));
}


#' Once the prior on the numbuer of mixture $M$ is assumed to be a Shifted Posson of parameter \code{Lambda}, 
#' this function  adopt a \emph{bisection method} to find the value of \code{gamma} such that the induced distribution
#'  on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use a bisection
#'   method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} 
#'   and an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.
#'     A defaault value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31),
#'      the function stops warning that convergence has not bee reached.
#'
#' @param n             The sample size
#' @param Lambda        The parameter of the Shifted Poisson for the number of components of the mixture
#' @param Kstar         The mean number of cluster the user want to specify
#' @param gam_min       The lower bound of the interval in which \code{gamma} should be lie
#' @param gam_max       The upper bound of the interval in which \code{gamma} should lie
#' @param tolerance     Tolerance of the method
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' Lam  <- 11
#' gam_po <-  AM_find_gamma_Pois(n,Lam,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1)
#' prior_K_po <-  AM_prior_K_Pois(n,gam_po,Lam)
#' prior_K_po\%*\%1:n

AM_find_gamma_Pois <- function (n,Lambda,Kstar=6, gam_min=0.0001,gam_max=10, tolerance=0.1) {
	return (find_gamma_Pois(n,Lambda,Kstar, gam_min,gam_max, tolerance));
}

#' Once the prior on the numbuer of mixture $M$ is assumed to be a Negative Binomial  Negative Binomial with parameter \code{r>0} and \code{0<p<1}, with  mean is 1+ r*p/(1-p), this function  adopt a \emph{bisection method} to find the value of \code{gamma} such that the induced distribution on the number of clusers is centered around a user specifed value \eqn{K^*}, i.e. the function use a bisection method to solve Eq.~{eq:findgamma} of WE NEED TO CITE ANTMAN PAPER. The user can provide a lower \eqn{\gamma_{l}} and an upper \eqn{\gamma_{u}} bound for the possible values of $gamma$. The default values are \eqn{\gamma_l= 10^{-3}} and \eqn{\gamma_{u}=10}.  A defaault value for the tolerance is \eqn{\epsilon=0.1}. Moreover, after a maximum number of iteration (default is 31), the function stops warning that convergence has not bee reached.
#'
#' @param n             The sample size
#' @param r      The dispersion parameter \code{r} of Negative Binomial
#' @param p      The probability of failure parameter \code{p} of Negative Binomial
#' @param Kstar         The mean number of cluster the user want to specify
#' @param gam_min  The lower bound of the interval in which \code{gamma} should be lie
#' @param gam_max   The upper bound of the interval in which \code{gamma} should lie
#' @param tolerance   tolerance of the method
#'
#'
#' @return A value of \code{gamma} such that \code{E(K)=K^*} 
#'
#' @keywords prior number of cluster
#'
#' @export
#' 
#' @examples
#' n <- 82
#' r <- 1
#' p <- 0.8571
#' gam_nb= AM_find_gamma_NegBin(n,r,p,Kstar=6, gam_min=0.001,gam_max=10000, tolerance=0.1)
#' prior_K_nb= AM_prior_K_NegBin(n,gam_nb, r, p)
#' prior_K_nb\%*\%1:n

AM_find_gamma_NegBin <- function (n,r,p,Kstar=6, gam_min=0.001,gam_max=10000, tolerance=0.1){
	return (find_gamma_NegBin(n,r,p,Kstar, gam_min,gam_max, tolerance));
}


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


#' this function compute the hyperparameters of an Normal-Inverse-Gamma 
#' distribution using an empirical Bayes approach.
#' @param  y       The data y
#' @param  scEmu a positive value scEmu (default =1) such that marginally E(mu)=(sample variance)*scEmu
#' @param  scEsig2 a positive value scEsig2 (default=3) such that marginally E(sig2)=(sample variance)*scEsig2
#' @param  CVsig2  The coefficient of variation of sig2 (default =3), CVsig2
#' @return Parameters of the normal inverse Gamma prior according to the empirical Bayes approach. 
#' @export
#' 

AM_emp_bayes_uninorm <- function(y,scEmu=1,scEsig2=3,CVsig2=3){
	n <- length(y)   ### sample size
	bary <- mean(y)  ### sample mean
	s2y <- var(y)    ### sample variance
	
	Emu <- bary
	Vmu <- s2y*scEmu
	Esig2 <- s2y/scEsig2
	Vsig2 <- CVsig2^2*Esig2^2
	
	return(from_statitstics2hyperparam(Emu,Vmu,Esig2,Vsig2))
	
}
