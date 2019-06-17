// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <limits>

#include <exception>      // to terminate a function by the command throw
// throw std::runtime_error("message")
//also includes std::exception, std::terminate


// [[Rcpp::export]]

Rcpp::NumericVector calcola_stirling(double gamma, int n){


	Rcpp::NumericVector out(n); // The output vector

	double app; //just for  my convenience in programming
	double lgammak; // again just for my convenience
	Rcpp::NumericVector Lgammamjg_over_Lgammajg(n); // again just for my convenience

	for(int k=1;k<=n;k++){
		Lgammamjg_over_Lgammajg[k-1]=std::lgamma(n+k*gamma)-std::lgamma(k*gamma);
		lgammak=std::lgammaf(k+1);

		out[k-1]=0;
		for(int j=1;j<=k;j++){
			app = R::lchoose(k,j)+Lgammamjg_over_Lgammajg[j-1]-lgammak;
			out[k-1] += std::pow(-1,j-k)*std::exp(app);
		}
	}

	return(out);

}


// [[Rcpp::export]]

Rcpp::NumericVector calcola_stirling_ricor(double gamma, unsigned int n){
	gamma=-gamma;
	Rcpp::NumericVector row_j(n+1,0.0); // The output vector initialize all the element to zero
	row_j[0]=1; /// Row j=0

	Rcpp::NumericVector row_jp1(n+1,0.0); // The output vector initialize all the element to zero
	row_jp1[1]=gamma; /// Row j+1=1



	for(unsigned int j=1;j<n;j++){

		std::copy(row_jp1.begin(),row_jp1.end(),row_j.begin()); // row j+1 becomes row j!
		//Rcpp::Rcout<<"j="<<j<<" row_j="<<row_j<<"\n";


		for(unsigned int k=1;k<=(j+1);k++){
			row_jp1[k]=(gamma*k-j)*row_j[k]+gamma*row_j[k-1];
		}

	}

	return(row_jp1[Rcpp::Range(1,n)]);

}



// [[Rcpp::export]]

Rcpp::NumericVector calcola_stirling_ricor_abs(double gamma, unsigned int n){

	Rcpp::NumericVector row_j(n+1,0.0); // The output vector initialize all the element to zero
	row_j[0]=1; /// Row j=0

	Rcpp::NumericVector row_jp1(n+1,0.0); // The output vector initialize all the element to zero

	row_jp1[0]=0; /// Row j+1=1
	row_jp1[1]=gamma; /// Row j+1=1



	for(unsigned int j=1;j<n;j++){

		std::copy(row_jp1.begin(),row_jp1.end(),row_j.begin()); // row j+1 becomes row j!
		//Rcpp::Rcout<<"j="<<j<<" row_j="<<row_j<<"\n";

		//row_jp1[0]=std::exp(std::lgamma(j-1+1)-std::lgamma(j-1-j+1));
		for(unsigned int k=1;k<=(j+1);k++){
			row_jp1[k]=(gamma*k+j)*row_j[k]+gamma*row_j[k-1];
		}

	}

	return(row_jp1[Rcpp::Range(1,n)]);

}








// [[Rcpp::export]]

Rcpp::NumericVector calcola_stirling_ricor_log(double gamma, unsigned int n){

	double infinito = std::numeric_limits<double>::infinity();

	Rcpp::NumericVector lrow_j(n+1,-infinito); // The output vector initialize all the element to zero
	lrow_j[0]=0; /// lRow j=0

	Rcpp::NumericVector lrow_jp1(n+1,-infinito); // The output vector initialize all the element to zero
	lrow_jp1[1]=std::log(gamma); /// lRow j+1=1



	for(unsigned int j=1;j<n;j++){

		std::copy(lrow_jp1.begin(),lrow_jp1.end(),lrow_j.begin()); // lrow j+1 becomes lrow j!
		//Rcpp::Rcout<<"j="<<j<<" lrow_j="<<lrow_j<<"\n";

		//Rcpp::Rcout<<"exp(-inf)="<< std::exp(lrow_j[0]-lrow_j[1])<<"\n";



		for(unsigned int k=1;k<=(j);k++){
			lrow_jp1[k]= std::log(gamma*k+j)+lrow_j[k]+ std::log(1+gamma/(gamma*k+j)*std::exp(lrow_j[k-1]-lrow_j[k]));
		}

		lrow_jp1[j+1]= (j+1)*std::log(gamma);

	}

	return(lrow_jp1[Rcpp::Range(1,n)]);

}












////////////////////////////////////////////////////////////////////////////
////////////////  Vnk FUNCTIONs     ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////











// [[Rcpp::export]]
Rcpp::NumericVector VnkPoisson(unsigned int n,double Lambda,double gamma){

	Rcpp::NumericVector out(n);
	int m;

	double ldenspoi;
	double controllo;
	for(int k=1;k<=n;k++){
		controllo=1.;
		m=0;
		out[k-1]=0;
		int lfactn=std::lgammaf(n);
		//Rcpp::Rcout<<"k="<<k<<"\n";

		Rcpp::NumericVector appoggio(0);
		while( (controllo>-1500) || (m<1500)  ){

			double first= std::lgammaf(m+k+1)-std::lgammaf(m+1);
			//Rcpp::Rcout<<"first= "<<first <<"\n";		

			ldenspoi=-Lambda+(m+k-1)*std::log(Lambda)-std::lgammaf(m+k);
			//Rcpp::Rcout<<"ldensita= "<<ldenspoi <<"\n";		

			double third=std::lgamma((k+m)*gamma) -std::lgamma((k+m)*gamma+n);
			//Rcpp::Rcout<<"third= "<<third <<"\n";		

			controllo=first+ldenspoi+third;
			//Rcpp::Rcout<<"m="<<m<<" controllo="<<controllo<<"\n";
			appoggio.push_back(controllo);
			m +=1;
		}
		double massimo=Rcpp::max(appoggio);
		Rcpp::NumericVector sottraggo= appoggio-massimo;
		//Rcpp::Rcout<<"sottraggo="<<Rcpp::sum(Rcpp::exp( sottraggo))<<"\n";

		out[k-1]=massimo+std::log(Rcpp::sum(Rcpp::exp( sottraggo)));
		//Rcpp::Rcout<<"out[k-1]="<<out[k-1]<<"\n";

	}
	return(out) ;
}



















// [[Rcpp::export]]
Rcpp::NumericVector VnkNegBin(unsigned int n,double r,double p,double gamma){

	Rcpp::NumericVector out(n);
	int m;

	double ldensnegbin;
	double controllo;
	for(unsigned int k=1;k<=n;k++){
		controllo=1.;
		m=0;
		out[k-1]=0;
		int lfactn=std::lgammaf(n);
		Rcpp::Rcout<<"k="<<k<<"\n";

		Rcpp::NumericVector appoggio(0);
		while( (controllo>-100) || (m<100)  ){

			double first= std::lgammaf(m+k+1)-std::lgammaf(m+1);
			//Rcpp::Rcout<<"first= "<<first <<"\n";		

			ldensnegbin=std::lgamma(r+m+k-1)-std::lgamma(r)-std::lgammaf(m+k)+(m+k-1)*std::log(p)+r*std::log(1-p);
			//Rcpp::Rcout<<"ldensita= "<<ldensgamma <<"\n";		

			double third=std::lgamma((k+m)*gamma) -std::lgamma((k+m)*gamma+n);
			//Rcpp::Rcout<<"third= "<<third <<"\n";		

			controllo=first+ldensnegbin+third;
			//Rcpp::Rcout<<"m="<<m<<" controllo="<<controllo<<"\n";
			appoggio.push_back(controllo);
			m +=1;
		}
		double massimo=Rcpp::max(appoggio);
		Rcpp::NumericVector sottraggo= appoggio-massimo;
		Rcpp::Rcout<<"sottraggo="<<Rcpp::sum(Rcpp::exp( sottraggo))<<"\n";

		out[k-1]=massimo+std::log(Rcpp::sum(Rcpp::exp( sottraggo)));
		Rcpp::Rcout<<"out[k-1]="<<out[k-1]<<"\n";

	}
	return(out) ;
}










// [[Rcpp::export]]
Rcpp::NumericVector VnkDelta(unsigned int n,unsigned int Mstar,double gamma){

	Rcpp::NumericVector out(n);
	int m;
	double infinito = std::numeric_limits<double>::infinity();


	for(unsigned int k=1;k<=n;k++){
		m=0;

		if(k<=Mstar){
			//int lfactn=std::lgammaf(n);
			Rcpp::Rcout<<"k="<<k<<"\n";
			Rcpp::NumericVector appoggio(0);
			out[k-1]=std::lgammaf(Mstar+1)-std::lgammaf(Mstar-k+1)+std::lgamma(gamma*Mstar)-std::lgamma(n+gamma*Mstar);
		}else{
			out[k-1]=-infinito;
		}
		Rcpp::Rcout<<"out[k-1]="<<out[k-1]<<"\n";

	}
	return(out) ;
}
















////////////////////////////////////////////////////////////////////////////
////////////////  PRIORS FOR K      ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////








// This function compute the prior on the number of cluster, i.e. occupied component of the mixutre
// under shifted Poisson for M.


// [[Rcpp::export]]
Rcpp::NumericVector prior_K_Pois(unsigned int n,double gamma,double Lambda){
	
	// Compute the log(V(m,k)) under the Poissom
	Rcpp::NumericVector vvv=VnkPoisson(n,Lambda,gamma);

	// Compute the Stirling number in log scale
	Rcpp::NumericVector stir=calcola_stirling_ricor_log(gamma, n);

	// Compute the prior for K
	Rcpp::NumericVector pstrk = Rcpp::exp(vvv+stir);
	
	double sum=Rcpp::sum(pstrk);
	if(std::abs(sum-1)>0.01){
		Rcpp::Rcout<<" Sorry I was unable to compute the prior on the number of cluster"<<"\n"<<"for the parameters"<<" n="<<n<<" gamma="<<gamma<<" Lambda="<<Lambda<<"\n";
		Rcpp::stop("Unexpected condition occurred");
	}
	
	for(unsigned int l=0;l<n;l++){
			pstrk[l]=pstrk[l]/sum;
	}
	return(pstrk);
}







// This function compute the prior on the number of cluster, i.e. occupied component of the mixutre
// under the Negative binomial for M.

// [[Rcpp::export]]
Rcpp::NumericVector prior_K_NegBin(unsigned int n,double gamma,double r, double p){
	
	// Compute the log(V(m,k)) under the Poissom
	Rcpp::NumericVector vvv= VnkNegBin(n,r,p,gamma);

	// Compute the Stirling number in log scale
	Rcpp::NumericVector stir=calcola_stirling_ricor_log(gamma, n);

	// Compute the prior for K
	Rcpp::NumericVector pstrk = Rcpp::exp(vvv+stir);
	
	double sum=Rcpp::sum(pstrk);
	if(std::abs(sum-1)>0.1){
		Rcpp::Rcout<<" Sorry I was unable to compute the prior on the number of cluster"<<"\n"<<"for the parameters"<<" n="<<n<<" gamma="<<gamma<<" r="<<r << " p="<<p<< "the mass is"<<sum<<"\n";
		Rcpp::stop("Unexpected condition occurred");
	}
	
	for(unsigned int l=0;l<n;l++){
			pstrk[l]=pstrk[l]/sum;
	}
	return(pstrk);
}






// This function compute the prior on the number of cluster, i.e. occupied component of the mixutre
// under the Delta prior for K.

// [[Rcpp::export]]
Rcpp::NumericVector prior_K_Delta(const unsigned int n,const double gamma,const unsigned Mstar){
	
	// Compute the log(V(m,k)) under the Delta
	Rcpp::NumericVector vvv= VnkDelta(n,Mstar,gamma);

	// Compute the Stirling number in log scale
	Rcpp::NumericVector stir=calcola_stirling_ricor_log(gamma, n);

	// Compute the prior for K
	Rcpp::NumericVector pstrk = Rcpp::exp(vvv+stir);
	
	double sum=Rcpp::sum(pstrk);
	if(std::abs(sum-1)>0.01){
		Rcpp::Rcout<<" Sorry I was unable to compute the prior on the number of cluster"<<"\n"<<"for the parameters"<<" n="<<n<<" Mstar="<<Mstar<<"\n";
		Rcpp::stop("Unexpected condition occurred");
	}
	
	for(unsigned int l=0;l<n;l++){
			pstrk[l]=pstrk[l]/sum;
	}
	return(pstrk);
}












////////////////////////////////////////////////////////////////////////////
////////////////  BISECTION METHODS ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////








/// Utility function
double calcola_media(const Rcpp::NumericVector prob,const unsigned n){
	
	double out=0;
	for(unsigned int i =0;i<n;i++){
		out += (i+1)*prob[i];
	}

	return(out);

}




// [[Rcpp::export]]
double find_gamma_Pois(const unsigned int n,const double Lambda,const unsigned int Kstar, double gam_min,double gam_max,const double tollerance,const unsigned int max_iter=30){

	
	//Rcpp::IntegerVector unoton=Rcpp:seq(1,n);
	
	Rcpp::NumericVector p_min=prior_K_Pois(n,gam_min,Lambda);
	double Kmin=calcola_media(p_min,n);
	Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

	Rcpp::NumericVector p_max=prior_K_Pois(n,gam_max,Lambda);
	double Kmax=calcola_media(p_max,n);
	Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

	if( (Kmin-Kstar)>0){
		Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";
		Rcpp::stop("Sorry (Kmin-Kstar)>0, you should try with a smaller value of gam_min");
	} 
	if( (Kmax-Kstar)<0){
		Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";
		Rcpp::stop("Sorry (Kmax-Kstar)<0, you should try with a larger value of gam_max");
	}

	double gam_mean;
	Rcpp::NumericVector p_mean;
	double Kmean; 
	double obs_tol=Kmax-Kmin;
	unsigned niter=0;

	while( (obs_tol>=tollerance) && (niter<max_iter) ){
		//Compute the center
		gam_mean=(gam_min+gam_max)/2;	
		p_mean=prior_K_Pois(n,gam_mean,Lambda);
		Kmean=calcola_media(p_mean,n);
		Rcpp::Rcout<<"K_mean="<<Kmean<<" gam_mean="<<gam_mean<<"\n";

		//If the center leads a number of cluster larger than the target
		if(Kmean-Kstar>0){
			gam_max=gam_mean;
			Kmax=Kmean;
			Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

		}
		else{//if the center leads a number of clusters samaller than the target
			gam_min=gam_mean;
			Kmin=Kmean;
			Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

		}
		obs_tol=Kmax-Kmin;
		niter+=1;
		
	}

	if(niter>=max_iter){
		Rcpp::warning("Not converged increase maximum number of iteration, max_iter");
	}
	return((gam_min+gam_max)/2);
}














///Metodo bisezione Negative binomial

// [[Rcpp::export]]
double find_gamma_NegBin(const unsigned int n,const double r, const double p,const unsigned int Kstar, double gam_min,double gam_max,const double tollerance,const unsigned int max_iter=30){

	
	//Rcpp::IntegerVector unoton=Rcpp:seq(1,n);
	
	Rcpp::NumericVector p_min=prior_K_NegBin(n,gam_min,r,p);
	double Kmin=calcola_media(p_min,n);
	Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

	Rcpp::NumericVector p_max=prior_K_NegBin(n,gam_max,r,p);
	double Kmax=calcola_media(p_max,n);
	Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

	if( (Kmin-Kstar)>0){
		Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";
		Rcpp::stop("Sorry (Kmin-Kstar)>0, you should try with a smaller value of gam_min");
	} 
	if( (Kmax-Kstar)<0){
		Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";
		Rcpp::stop("Sorry (Kmax-Kstar)<0, you should try with a larger value of gam_max");
	}

	double gam_mean;
	Rcpp::NumericVector p_mean;
	double Kmean; 
	double obs_tol=Kmax-Kmin;
	unsigned niter=0;

	while( (obs_tol>=tollerance) && (niter<max_iter) ){
		//Compute the center
		gam_mean=(gam_min+gam_max)/2;	
		p_mean=prior_K_NegBin(n,gam_mean,r,p);
		Kmean=calcola_media(p_mean,n);
		Rcpp::Rcout<<"K_mean="<<Kmean<<" gam_mean="<<gam_mean<<"\n";

		//If the center leads a number of cluster larger than the target
		if(Kmean-Kstar>0){
			gam_max=gam_mean;
			Kmax=Kmean;
			Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

		}
		else{//if the center leads a number of clusters samaller than the target
			gam_min=gam_mean;
			Kmin=Kmean;
			Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

		}
		obs_tol=Kmax-Kmin;
		niter+=1;
		
	}

	if(niter>=max_iter){
		Rcpp::warning("Not converged increase maximum number of iteration, max_iter");
	}
	return((gam_min+gam_max)/2);
}










///Metodo bisezione Delta

// [[Rcpp::export]]
double find_gamma_Delta(const unsigned int n,const unsigned Mstar,const unsigned int Kstar, double gam_min,double gam_max,const double tollerance,const unsigned int max_iter=30){

	//TODO: DFai il check di Kstar minore di M_star e n
	
	//Rcpp::IntegerVector unoton=Rcpp:seq(1,n);
	
	Rcpp::NumericVector p_min= prior_K_Delta(n,gam_min,Mstar);
	double Kmin=calcola_media(p_min,n);
	Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

	Rcpp::NumericVector p_max=prior_K_Delta(n,gam_max,Mstar);
	double Kmax=calcola_media(p_max,n);
	Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

	if( (Kmin-Kstar)>0){
		Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";
		Rcpp::stop("Sorry (Kmin-Kstar)>0, you should try with a smaller value of gam_min");
	} 
	if( (Kmax-Kstar)<0){
		Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";
		Rcpp::stop("Sorry (Kmax-Kstar)<0, you should try with a larger value of gam_max");
	}

	double gam_mean;
	Rcpp::NumericVector p_mean;
	double Kmean; 
	double obs_tol=Kmax-Kmin;
	unsigned niter=0;

	while( (obs_tol>=tollerance) && (niter<max_iter) ){
		//Compute the center
		gam_mean=(gam_min+gam_max)/2;	
		p_mean=prior_K_Delta(n,gam_mean,Mstar);
		Kmean=calcola_media(p_mean,n);
		Rcpp::Rcout<<"K_mean="<<Kmean<<" gam_mean="<<gam_mean<<"\n";

		//If the center leads a number of cluster larger than the target
		if(Kmean-Kstar>0){
			gam_max=gam_mean;
			Kmax=Kmean;
			Rcpp::Rcout<<"K_max="<<Kmax<<" gam_max="<<gam_max<<"\n";

		}
		else{//if the center leads a number of clusters samaller than the target
			gam_min=gam_mean;
			Kmin=Kmean;
			Rcpp::Rcout<<"K_min="<<Kmin<<" gam_min="<<gam_min<<"\n";

		}
		obs_tol=Kmax-Kmin;
		niter+=1;
		
	}

	if(niter>=max_iter){
		Rcpp::warning("Not converged increase maximum number of iteration, max_iter");
	}
	return((gam_min+gam_max)/2);
}









//////////////////////////////////////////////////////////////////////////////////////////


