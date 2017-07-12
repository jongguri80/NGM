
// includes from the plugin
#include <RcppEigen.h>
#include <Rcpp.h>
#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif
#ifndef END_RCPP
#define END_RCPP
#endif
using namespace Rcpp;
// user includes
// declarations
extern "C"{
SEXP rcpp_loglik( SEXP Y_obs, SEXP del_obs, SEXP V_eps, SEXP V_xi, SEXP Rho, SEXP mu_int, SEXP C_int) ;
}
// definition
SEXP rcpp_loglik( SEXP Y_obs, SEXP del_obs, SEXP V_eps, SEXP V_xi, SEXP Rho, SEXP mu_int, SEXP C_int ){
	BEGIN_RCPP
	using Eigen::Map;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Rcpp::as;
	typedef Map<VectorXd> MapVecd;
	typedef Map<MatrixXd> MapMatd;

	const MapVecd Y(as<MapVecd > (Y_obs));
	const MapVecd del(as<MapVecd> (del_obs));
	
	const double sigma2_eps = as<double>(V_eps);
	const double sigma2_xi = as<double>(V_xi);
	const double rho = as<double>(Rho);  
	
	const MapVecd mu00 (as<MapVecd> (mu_int));
	const MapMatd C00 (as<MapMatd> (C_int));

	int J(Y.size()); 

	VectorXd a, mu0;
	MatrixXd FF, GG, C0, W, R, A; 
	mu0 = mu00; C0 = C00;
	
	GG = MatrixXd::Identity(3, 3); W = MatrixXd::Zero(3, 3); 
	FF = MatrixXd::Zero(1, 3);FF(0,0) = 1.0;
	double x1, x2, Q, e, f, loglik;
	loglik = 0.0;
	for(int j=0; j < J; ++j){
		x1 = exp(-rho*del(j)); x2 = exp(-2.0*rho*del(j));
		GG(0,1) = (1.0-x1)/rho; GG(0,2) = del(j) - GG(0,1);
		GG(1,1) = x1; GG(1,2) = 1.0 - x1;
			
		a = GG*mu0;  f = (FF*a).sum();  e = Y(j) - f; 
		
		W(0,0) = del(j)/(rho*rho) + 0.5*(-3.0 + 4.0*x1 - x2)/(rho*rho*rho);
		W(0,1) = W(1,0) = 0.5*(1.0 - 2.0*x1 + x2)/(rho*rho);
		W(1,1) = 0.5*(1.0-x2)/rho;		

		R = GG*C0*GG.transpose() + sigma2_xi*W;
		Q = (FF*R*FF.transpose()).sum() + sigma2_eps;

		loglik += 0.5*(log(Q) + e*e/Q);

		// update priors //
		A = R*FF.transpose()/Q;
		mu0 = a + A*e;
		C0 = R - A*A.transpose()*Q;
	}
	return Rcpp::wrap(loglik);

	END_RCPP
}
