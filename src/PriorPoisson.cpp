/*
 * PoissonPrior.cpp
 *
 *  Created on: Apr 12, 2019
 *      Author: toky
 */


#include "PriorPoisson.hpp"

#include <cmath>
#include <vector>



void PriorPoisson::update_q_param (const  double U, const  int K) {
		double lphi_u=-h_param.gamma*std::log(1+U);
		double lpeso=lphi_u-std::log(1+bq);
		double lunif=std::log(R::runif(0.0,1.1));
		double astar=K+aq ;// See the notation of Point 4 in 10.1
		double rate=(1-std::exp(lphi_u)+bq);
		this->q_param.lambda = lunif<lpeso ? R::rgamma(astar+1,1.0/rate) : R::rgamma(astar,1.0/rate);
}

inline double log_full_gamma( const double Loc_gamma, const int K ,const  std::vector<int> & nj,const   double Lambda_current,const  double U_current ,const  double ag,const  double bg){

	double out=0;
	double up1g=std::pow(1+U_current,Loc_gamma);

	out+=std::log(Lambda_current/up1g+K)+Lambda_current/up1g-K*std::log(up1g);
	for(int j=0;j<K;j++){
		out+=std::lgamma(Loc_gamma+ (double) nj[j])-std::lgamma(Loc_gamma);
	}
	/// When the prior is a gamma

	out+=(ag-1)*std::log(Loc_gamma)-bg*Loc_gamma;

	return(out);
}

void PriorPoisson::update_h_param (const  double U, const  int K, const std::vector<int> &nj ) {

		const double vecchio = this->h_param.gamma;
		const double lmedia = std::log(vecchio);

		//Propose a new value
		const double lnuovo=R::rnorm(lmedia,lsd);
		const double nuovo=std::exp(lnuovo);

		double ln_acp = log_full_gamma(nuovo , K , nj,  this->q_param.lambda, U , ah, bh ) - lmedia;

		ln_acp= ln_acp - (log_full_gamma(vecchio , K , nj,   this->q_param.lambda, U , ah, bh) - lnuovo);

		const double lnu=std::log(R::runif(0.0,1.0));


		this->h_param.gamma = lnu<ln_acp ? nuovo : vecchio;

}


void PriorPoisson::update_h_param_adapitive (const int iter, const  double U, const  int K, const std::vector<int> &nj ) {

		const double vecchio = this->h_param.gamma;
		const double lmedia = std::log(vecchio);

		//Propose a new value
		const double lnuovo=R::rnorm(lmedia,lsd);
		const double nuovo=std::exp(lnuovo);

		double ln_acp = log_full_gamma(nuovo , K , nj,  this->q_param.lambda, U , ah, bh ) - lmedia;

		ln_acp= ln_acp - (log_full_gamma(vecchio , K , nj,   this->q_param.lambda, U , ah, bh) - lnuovo);

		const double lnu=std::log(R::runif(0.0,1.0));


		this->h_param.gamma = lnu<ln_acp ? nuovo : vecchio;


		// This is a new parameter to adjust lsd (ADAPTIVE METROPOLIS; the user could be allowed to set a different value in (-1,0) different that -0.7; Even if it is dangerous to change it 
		double wg=std::pow(iter+1,-0.7);
		// This is a new parameter to adjust lsd (ADAPTIVE METROPOLIS; the user could be allowed to set a different value in (0,1) different that 0.234; We should worn however. 
		// The adaptive Rejection Metropolis Hasting we are going to use is Algorithm 5 of Griffin Stephens (2003)
		double bartau =0.234;

		// Adaptive metropolis : Algorithm 5 Griffin Sthephens 
		//will this change lsd for the next iteration?
		lsd= lsd+wg*(std::exp(std::min(0.0,ln_acp))-bartau);
		if(lsd<std::pow(10,-50)){
			lsd=std::pow(10,-50);
		}
		if(lsd>std::pow(10,50)){
			lsd=std::pow(10,50);
		}
		std::cout<<"lsd= "<<lsd<<"\n";


}








int PriorPoisson::init_M_na() {
	return R::rpois(this->q_param.lambda);
}
int PriorPoisson::update_M_na(const double U ,  const int K) {

	int M_na;

	const double phi_u=-this->h_param.gamma*log(1+U);
	const double Lambda_u=exp(std::log(this->q_param.lambda)+phi_u);

	const double unif=R::runif(0.0,1.0);

	if(unif<(Lambda_u/(Lambda_u+K))){
		M_na=R::rpois(Lambda_u)+1;

	}
	else{
		M_na=R::rpois(Lambda_u);
	}
	return M_na;

}
