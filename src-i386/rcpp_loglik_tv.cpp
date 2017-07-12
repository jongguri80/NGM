
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
SEXP rcpp_loglik_tv( SEXP Y_obs, SEXP t_obs, SEXP V_eps, SEXP Rho, SEXP Nu0, SEXP Nu1, SEXP V_xi, SEXP Mu00, SEXP C00, SEXP T0) ;
}
// definition
SEXP rcpp_loglik_tv( SEXP Y_obs, SEXP t_obs, SEXP V_eps, SEXP Rho, SEXP Nu0, SEXP Nu1, SEXP V_xi, SEXP Mu00, SEXP C00, SEXP T0){
	BEGIN_RCPP
	using Eigen::Map;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Rcpp::as;
	typedef Map<VectorXd> MapVecd;
	typedef Map<MatrixXd> MapMatd;

	const MapVecd Y(as<MapVecd > (Y_obs));
	const MapVecd tt(as<MapVecd> (t_obs));
	
	const double sigma2_eps = as<double>(V_eps);
	const double rho = as<double>(Rho);  
	const double nu0 = as<double>(Nu0);
	const double nu1 = as<double>(Nu1);
	const double sigma2_xi = as<double>(V_xi);
	const double t0 = as<double>(T0);
	
	const MapVecd mu00 (as<MapVecd> (Mu00));
	const MapMatd c00 (as<MapMatd> (C00));

	int J(Y.size()); 
	
	VectorXd a, mu0, del; 
	mu0 = VectorXd::Zero(4); del = VectorXd::Zero(J);
	mu0.segment(0,2) = mu00; mu0(2) = nu0; mu0(3) = nu1;
	del(0) = tt(0)-t0; for(int j=1;j<J;++j) del(j) = tt(j) - tt(j-1);

	MatrixXd FF, GG, C0, W, R, A; 
    
    C0 = MatrixXd::Zero(4,4); C0.block(0,0,2,2) = c00;
	
	GG = MatrixXd::Zero(2, 4); W = MatrixXd::Zero(2, 2); 
	FF = MatrixXd::Zero(1, 2);FF(0,0) = 1.0;
	double x1, x2, x3, Q, e, f, loglik, t_old;
	loglik = 0.0;
	
	for(int j=0; j < J; ++j){
		x1 = exp(-rho*del(j)); x2 = (1.0-x1)/rho; x3 = exp(-2.0*rho*del(j));
		t_old = tt(j) - del(j);
			
		GG(0,0) = 1.0; 
        GG(0,1) = x2; 
        GG(0,2) = del(j) - x2; 
        GG(0,3) = .5*del(j)*del(j) + t_old*GG(0,2);
        
		GG(1,1) = x1; 
        GG(1,2) = (1.0 - x1); 
        GG(1,3) = del(j) + t_old*GG(1,2);
			
		a = GG*mu0; // 2x4, 4x1 -> 2x1
        f = (FF*a).sum();  // 1x2, 2x1 -> d
        e = Y(j) - f; //d
		
		W(0,0) = del(j)/(rho*rho) + 0.5*(-3.0 + 4.0*x1 - x3)/(rho*rho*rho);
		W(0,1) = W(1,0) = 0.5*(1.0 - 2.0*x1 + x3)/(rho*rho);
		W(1,1) = 0.5*(1.0-x3)/rho;		

		R = GG*C0*GG.transpose() + sigma2_xi*W; // 2x4, 4x4, 4x2 -> dim(R) = 2x2
		Q = (FF*R*FF.transpose()).sum() + sigma2_eps; // 1x2, 2x2, 2x1 -> dim(Q) = 1x1

		loglik += 0.5*(log(Q) + e*e/Q);

		// update priors //
		A = R*FF.transpose()/Q; // 2x2, 2x1 -> 2x1
		mu0.segment(0,2) = a + A*e; 
		C0.block(0,0,2,2) = R - A*A.transpose()*Q;
	}
	return Rcpp::wrap(loglik);

	END_RCPP
}





