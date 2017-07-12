
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
SEXP rcpp_sample_mu( SEXP Y_obs, SEXP del_obs, SEXP V_eps, SEXP V_xi, SEXP Rho, SEXP mu_int, SEXP C_int, SEXP Z_mat) ;
}
// definition
SEXP rcpp_sample_mu( SEXP Y_obs, SEXP del_obs, SEXP V_eps, SEXP V_xi, SEXP Rho, SEXP mu_int, SEXP C_int, SEXP Z_mat){
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
	const MapMatd Z(as<MapMatd> (Z_mat));
	
	
	int J(Y.size()); 

	VectorXd a, mu0, h;
	MatrixXd FF, GG, C0, W, iW, V0, tmp_mat;
	mu0 = mu00; C0 = C00;
	
	GG = MatrixXd::Identity(3, 3); iW = W = MatrixXd::Zero(3, 3); 
	FF = MatrixXd::Zero(1, 3);FF(0,0) = 1.0;
	double x1, x2, e, f;

	MatrixXd mat_a, mat_m; 
	MatrixXd* mat_iC0; mat_iC0 = new MatrixXd[J+1];
	
	mat_a = MatrixXd::Zero(J, 2);  mat_m = MatrixXd::Zero(J+1, 2);
	mat_m.row(0) = mu0.segment(0,2);
	mat_iC0[0] = C0.block(0,0,2,2).inverse();
	
	for(int j=0; j < J; ++j){
		x1 = exp(-rho*del(j)); x2 = exp(-2.0*rho*del(j));
		GG(0,1) = (1.0-x1)/rho; GG(0,2) = del(j) - GG(0,1);
		GG(1,1) = x1; GG(1,2) = 1.0 - x1;
			
		a = GG*mu0; 
		f = (FF*a).sum(); 
		e = Y(j) - f; 
		
		W(0,0) = del(j)/(rho*rho) + 0.5*(-3.0 + 4.0*x1 - x2)/(rho*rho*rho);
		W(0,1) = W(1,0) = 0.5*(1.0 - 2.0*x1 + x2)/(rho*rho);
		W(1,1) = 0.5*(1.0-x2)/rho;		
		
		V0 = GG*C0*GG.transpose() + sigma2_xi*W;

		// update priors //
		tmp_mat = V0.block(0,0,2,2).inverse() + FF.block(0,0,1,2).transpose()*FF.block(0,0,1,2)/sigma2_eps;
		C0.block(0,0,2,2) = tmp_mat.inverse();
		mu0 = a + C0*(FF.transpose()*e/sigma2_eps); 
		// save parameters //
		mat_a.row(j) = a.segment(0,2);
		mat_m.row(j+1) = mu0.segment(0,2);
		mat_iC0[j+1] = tmp_mat;
	}
	
	MatrixXd Vmat; Vmat = MatrixXd::Zero(2,2);
	MatrixXd mat_theta, mat_h; mat_theta = MatrixXd::Zero(J+1, 2);
	mat_h = MatrixXd::Zero(J, 2);
	
	MatrixXd* mat_Vmat; mat_Vmat = new MatrixXd[J];
	VectorXd diff;
	Eigen::LLT<MatrixXd> llt(C0.block(0,0,2,2));
	mat_theta.row(J) = mat_m.row(J) + (MatrixXd(llt.matrixL())*Z.row(J).transpose()).transpose(); 
	
	for(int j=J-1; j >= 0; --j){
		x1 = exp(-rho*del(j)); x2 = exp(-2.0*rho*del(j));
		GG(0,1) = (1.0-x1)/rho; GG(0,2) = del(j) - GG(0,1);
		GG(1,1) = x1; GG(1,2) = 1.0 - x1;
			
		a = GG*mu0; f = (FF*a).sum(); e = Y(j) - f; 
		
		W(0,0) = del(j)/(rho*rho) + 0.5*(-3.0 + 4.0*x1 - x2)/(rho*rho*rho);
		W(0,1) = W(1,0) = 0.5*(1.0-2.0*x1 + x2)/(rho*rho);
		W(1,1) = 0.5*(1.0-x2)/rho;	
		W = sigma2_xi*W;
		iW.block(0,0,2,2) = W.block(0,0,2,2).inverse();
		

		tmp_mat = mat_iC0[j] + GG.block(0,0,2,2).transpose()*iW.block(0,0,2,2)*GG.block(0,0,2,2);
		Vmat = tmp_mat.inverse();
		diff = mat_theta.row(j+1) - mat_a.row(j);
		
		h = mat_m.row(j) +  (Vmat*(GG.block(0,0,2,2).transpose()*iW.block(0,0,2,2))*diff).transpose();
		
		mat_Vmat[j] = Vmat;
		mat_h.row(j) = h;
		Eigen::LLT<MatrixXd> llt2(Vmat);
		mat_theta.row(j) = h + (MatrixXd(llt2.matrixL())*Z.row(j).transpose());
	}	
	
	delete [] mat_Vmat;
	delete [] mat_iC0;
	return Rcpp::wrap(mat_theta);	

	END_RCPP
}
