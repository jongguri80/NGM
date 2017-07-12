
SSVMOU_tv_Y_Pred <- function(theta_old, t_old, del, sigma2_eps, rho, nu0, nu1, sigma2_xi){
		x1 = exp(-rho*del); x2 = (1-x1)/rho; x3 = exp(-2*rho*del);
		GG = matrix(0, nrow=2, ncol=4)
		GG[1,1] = 1
		GG[1,2] = x2
		GG[1,3] = del - x2
		GG[1,4] = .5*del^2 + t_old*GG[1,3]
		
		GG[2,2] = x1
		GG[2,3] = (1 - x1)
		GG[2,4] = del + t_old*GG[2,3]

		m = GG %*% theta_old
		
		W = diag(2)
		W[1,1] = del/rho^2 + 0.5*(-3 + 4*x1 - x3)/rho^3
		W[1,2] = W[2,1] = 0.5*(1-2*x1 + x3)/rho^2
		W[2,2] = 0.5*(1-x3)/rho
		W = sigma2_xi*W
			
		theta_new = c(t(chol(W)) %*% rnorm(2) + m, nu0, nu1)
		return(list(theta=theta_new, Y = rnorm(1, theta_new[1], sqrt(sigma2_eps)), mu = m, Sigma = W))
}


simul_SSVMtv = function(n_i, theta0, time_list, sigma2_eps, rhos, nu0, nu1, sigma2_xi, t0 = -0.01){
	N = length(n_i)
	theta0 = cbind(theta0, nu0, nu1)
	Y = theta = as.list(NULL)
	for(i in 1:N){
		del = diff(c(t0, time_list[[i]]))
		
		theta[[i]] = matrix(0, nrow = n_i[i], ncol = 4)
		Y[[i]] = rep(0, n_i[i])
		t00 = time_list[[i]][1] - del[1]
		tmp=SSVMOU_tv_Y_Pred(theta0[i,], t00, del[1], sigma2_eps, rhos[i], nu0[i], nu1[i], sigma2_xi)
		theta[[i]][1,] = tmp$theta
		Y[[i]][1] = tmp$Y
		for(j in 2:n_i[i]){
			tmp=SSVMOU_tv_Y_Pred(theta[[i]][j-1,], time_list[[i]][j-1], del[j], sigma2_eps, rhos[i], nu0[i], nu1[i], sigma2_xi)
			Y[[i]][j] = tmp$Y
			theta[[i]][j,] = tmp$theta
		}	
	}
	return(list(Y = Y, theta=theta))
}


loglik_tv_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=-0.01){
	y_i = Y_list[[i]]
	rho_i = rho[i]
	nu0_i = nu0[i]
	nu1_i = nu1[i]
	del_list = diff(c(t0, time_list[[i]]))
	
	mu0 = c(mu0[[i]], nu0_i, nu1_i)
	C0=diag(c(diag(C0[[i]]), 0, 0))

	GG = matrix(0, nrow = 2, ncol=4)
	FF = matrix(c(1, 0), 1, 2)
	loglik = 0
	W = matrix(0, 2, 2)
	for(j in 1:length(y_i)){
		del = del_list[j]
		t_old = time_list[[i]][j]-del_list[j]
		
		x1 = exp(-rho_i*del); x2 = (1-x1)/rho_i; x3 = exp(-2*rho_i*del);
		GG[1,1] = 1;
		GG[1,2] = x2;
		GG[1,3] = del - x2; 
		GG[1,4] = .5*del^2 + t_old*GG[1,3]
		
		GG[2,2] = x1;
		GG[2,3] = (1 - x1); 
		GG[2,4] = del + t_old*GG[2,3]	
		
		a = GG %*% mu0; f = FF %*% a; e = y_i[j] - f
		
		W[1,1] = del/rho_i^2 + 0.5*(-3 + 4*x1 - x3)/rho_i^3
		W[1,2] = W[2,1] = 0.5*(1-2*x1 + x3)/rho_i^2
		W[2,2] = 0.5*(1-x3)/rho_i
		
		V0 = GG %*% C0 %*% t(GG) + sigma2_xi*W
		Sigma = FF %*% V0 %*% t(FF) + sigma2_eps
		loglik = loglik + 0.5*(log(Sigma) + e^2/Sigma)
		
		## update priors ##
		C0[1:2, 1:2] = solve(solve(V0) + matrix(c(1/sigma2_eps, rep(0,3)), 2,2))
		mu0 = a + C0[1:2, 1:2] %*% (t(FF) %*% e/sigma2_eps)
		mu0 = c(mu0, nu0_i, nu1_i)
	}
	return(loglik)
}


Uni_loglik_tv2 = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, t0=-.01){
				N = length(Y_list)
				return(sum(vapply(1:N, loglik_tv_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, t0, FUN.VALUE = 0)))
}

cpp_loglik_tv_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=-.01){
				loglik = rcpp_loglik_tv(Y_obs=Y_list[[i]], t_obs=time_list[[i]], 
									V_eps = sigma2_eps, Rho=rho[i], Nu0 = nu0[i], Nu1 = nu1[i], V_xi = sigma2_xi, 
									Mu00=mu0[[i]], C00=C0[[i]], T0 = t0)
				return(loglik)
}


Uni_loglik_tv = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=-.01){
				N = length(Y_list)
				return(sum(vapply(1:N, cpp_loglik_tv_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0, FUN.VALUE = 0)))
}

############### sampling mean and velocity ########################3

cpp_mu_tv_i = function(i, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=-.01)
{
	J = length(time_list[[i]])
	Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
	mu = rcpp_sample_tv_mu(Y_obs = Y[[i]], t_obs = time_list[[i]], 
						   V_eps=sigma2_eps, Rho=rho[i], Nu0=nu0[i], Nu1= nu1[i], V_xi=sigma2_xi, 
						   Mu00=mu0[[i]], C00=C0[[i]], T0 = t0, Z_mat = Z)
	return(mu)
}


SampleMu_tv_cpp = function(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=-.01){
			  N = length(Y)
			  mu = lapply(X=1:N, cpp_mu_tv_i, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0, t0=t0)
			  return(mu)
}



#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
SS_tv_sigma2_eps = function(neg_loglik, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
						bounds=list(lower=0, upper=50), prior=list(a0=0.001, b0=0.001))
{
	
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_sigma2_eps = rinvgamma_trunc(1, a=lower, b=upper, shape=prior$a0, rate=prior$b0)
		neg_loglik = Uni_loglik_tv(Y_list, time_list, cand_sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		
		if(v > neg_loglik){sampled = 1; sigma2_eps = cand_sigma2_eps;}
		if( sigma2_eps < cand_sigma2_eps){upper = cand_sigma2_eps}else{lower = cand_sigma2_eps}
	}
	return(list(sigma2_eps=sigma2_eps, loglik=neg_loglik))
}


#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
SS_tv_sigma2_xi = function(neg_loglik, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
						bounds=list(lower=0, upper=50), prior=list(a0=.1, b0=.1))
{
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_sigma2_xi = rinvgamma_trunc(1, a=lower, b=upper, shape=prior$a0, rate=prior$b0)
		neg_loglik = Uni_loglik_tv(Y_list, time_list, sigma2_eps, rho, nu0, nu1, cand_sigma2_xi, mu0, C0)
		
		if(v > neg_loglik){sampled = 1; sigma2_xi = cand_sigma2_xi;}
		if( sigma2_xi < cand_sigma2_xi){upper = cand_sigma2_xi}else{lower = cand_sigma2_xi}
	}
	return(list(sigma2_xi=sigma2_xi, loglik=neg_loglik))
}


#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
SS_tv_nu1 = function(neg_loglik, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
						bounds=list(lower=-100, upper=100), prior=list(m0=0, v0=10^8))
{
	N = length(Y_list)
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_nu1 = rep(rtruncnorm(1, a=lower, b=upper, mean=prior$m0, sd=sqrt(prior$v0)),N)
		neg_loglik = Uni_loglik_tv(Y_list, time_list, sigma2_eps, rho, nu0, cand_nu1, sigma2_xi, mu0, C0)
		
		if(v > neg_loglik){sampled = 1; nu1 = cand_nu1;}
		if( nu1[1] < cand_nu1[1]){upper = cand_nu1[1]}else{lower = cand_nu1[1]}
	}
	return(list(nu1=nu1, loglik=neg_loglik))
}

SS_tv_nu0_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mu0, C0)
{

	neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	v = rexp(1) + neg_loglik_i 
	cand_nu0 = nu0
	sampled = 0; 
	bounds = mean_nu0[i] + 10*c(-1,1)*sqrt(sigma2_nu0);
	while(sampled == 0){
		cand_nu0[i] = rtruncnorm(1, a=bounds[1], b=bounds[2], mean=mean_nu0[i], sd=sqrt(sigma2_nu0))
		neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, cand_nu0, nu1, sigma2_xi, mu0, C0)

		if(v > neg_loglik_i){sampled = 1; nu0[i] = cand_nu0[i];}
		if( nu0[i] < cand_nu0[i]){bounds[2] = cand_nu0[i]}else{bounds[1] = cand_nu0[i]}
	}
	return(nu0[i])
}

SS_tv_nu0_par = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mu0, C0)
{
	N = length(Y_list)
	posterior = vapply(1:N, SS_tv_nu0_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mu0, C0, FUN.VALUE = 0)
	return(posterior)
}


SS_tv_nu1_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu1, sigma2_nu1, mu0, C0)
{
	neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	v = rexp(1) + neg_loglik_i 
	cand_nu1 = nu1
	sampled = 0; 
	bounds = mean_nu1[i] + 10*c(-1,1)*sqrt(sigma2_nu1);
	while(sampled == 0){
		cand_nu1[i] = rtruncnorm(1, a=bounds[1], b=bounds[2], mean=mean_nu1[i], sd=sqrt(sigma2_nu1))
		neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, cand_nu1, sigma2_xi, mu0, C0)

		if(v > neg_loglik_i){sampled = 1; nu1[i] = cand_nu1[i];}
		if( nu1[i] < cand_nu1[i]){bounds[2] = cand_nu1[i]}else{bounds[1] = cand_nu1[i]}
	}
	return(nu1[i])
}

SS_tv_nu1_par = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu1, sigma2_nu1, mu0, C0)
{
	N = length(Y_list)
	posterior = vapply(1:N, SS_tv_nu1_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu1, sigma2_nu1, mu0, C0, FUN.VALUE = 0)
	return(posterior)
}



SS_tv_phi_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, phi, sigma2_xi, mean_phi, sigma2_phi, mu0, C0)
{
	nu1 = nu0*phi
	neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	v = rexp(1) + neg_loglik_i 
	cand_phi = phi
	sampled = 0; 
	bounds = mean_phi[i] + 10*c(-1,1)*sqrt(sigma2_phi);
	while(sampled == 0){
		cand_phi[i] = rtruncnorm(1, a=bounds[1], b=bounds[2], mean=mean_phi[i], sd=sqrt(sigma2_phi))
		nu1 = nu0*cand_phi;
		neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)

		if(v > neg_loglik_i){sampled = 1; phi[i] = cand_phi[i];}
		if( phi[i] < cand_phi[i]){bounds[2] = cand_phi[i]}else{bounds[1] = cand_phi[i]}
	}
	return(phi[i])
}

SS_tv_phi_par = function(Y_list, time_list, sigma2_eps, rho, nu0, phi, sigma2_xi, mean_phi, sigma2_phi, mu0, C0)
{
	N = length(Y_list)
	posterior = vapply(1:N, SS_tv_phi_i, Y_list, time_list, sigma2_eps, rho, nu0, phi, sigma2_xi, mean_phi, sigma2_phi, mu0, C0, FUN.VALUE = 0)
	return(posterior)
}



SS_tv_nus_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mean_nu1, sigma2_nu1, mu0, C0)
{
	neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	v = rexp(1) + neg_loglik_i 
	cand_nu0 = nu0; cand_nu1 = nu1
	sampled = 0; 
	bounds0 = mean_nu0[i] + 10*c(-1,1)*sqrt(sigma2_nu0);
	bounds1 = mean_nu1[i] + 10*c(-1,1)*sqrt(sigma2_nu1);
	while(sampled == 0){
		cand_nu0[i] = rtruncnorm(1, a=bounds0[1], b=bounds0[2], mean=mean_nu0[i], sd=sqrt(sigma2_nu0))
		cand_nu1[i] = rtruncnorm(1, a=bounds1[1], b=bounds1[2], mean=mean_nu1[i], sd=sqrt(sigma2_nu1))

		neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, cand_nu0, cand_nu1, sigma2_xi, mu0, C0)

		if(v > neg_loglik_i){sampled = 1; nu0[i] = cand_nu0[i]; nu1[i] = cand_nu1[i];}
		if( nu0[i] < cand_nu0[i]){bounds0[2] = cand_nu0[i]}else{bounds0[1] = cand_nu0[i]}
		if( nu1[i] < cand_nu1[i]){bounds1[2] = cand_nu1[i]}else{bounds1[1] = cand_nu1[i]}
	}
	return(c(nu0[i], nu1[i]))
}

SS_tv_nus_par = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mean_nu1, sigma2_nu1, mu0, C0)
{
	N = length(Y_list)
	posterior = sapply(1:N, SS_tv_nus_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mean_nu1, sigma2_nu1, mu0, C0)
	return(t(posterior))
}


############## sample for rho ##############
SS_tv_rho_i = function(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0)
{
	neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	v = rexp(1) + neg_loglik_i
	cand_rho = rho
	sampled = 0;
	bounds = mean_logrho[i] + 10*c(-1,1)*sqrt(sigma2_logrho);
	while(sampled == 0){
		cand_rho[i] = exp(rtruncnorm(1, a=bounds[1], b=bounds[2], mean=mean_logrho[i], sd=sqrt(sigma2_logrho)))
		neg_loglik_i = cpp_loglik_tv_i(i, Y_list, time_list, sigma2_eps, cand_rho, nu0, nu1, sigma2_xi, mu0, C0)

		if(v > neg_loglik_i){sampled = 1; rho[i] = cand_rho[i]; }
		if( log(rho[i]) < log(cand_rho[i])){bounds[2] = log(cand_rho[i])}else{bounds[1] = log(cand_rho[i])}
	}
	return(rho[i])
}

SS_tv_rho_par = function(Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0)  
{
	N = length(Y_list)
	rhos = vapply(1:N, SS_tv_rho_i, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0, FUN.VALUE = 0)
	return(rhos)
}


######################################################    SSVM-OU ############################################################
BSSVMOU_tv = function(Y, X_logrho, X_nu0, time_list, mu0 = rep(list(c(0,0)), length(Y)), C0 = rep(list(diag(10^5, 2)),length(Y)), 
					smooth_fix = FALSE, lambda = 1,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.1, rho = rep(5, length(Y)), 
									nu0 = rep(0, length(Y)), nu1 = rep(0, length(Y)), 
									sigma2_logrho = .1, sigma2_nu0 = .1),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 200,
									  a0 = .1, b0 = .1,
									  betas0 = rep(0, ncol(X_logrho)), v_betas0 = rep(10^8, ncol(X_logrho)),
									  gammas0 = rep(0, ncol(X_nu0)), v_gammas0 = rep(10^8, ncol(X_nu0))),
					n_iter = 1500, burn_in=500, thin=1, per_print=500){
	#### assign initial parameters ####
	N = length(Y)
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	mean_nu0 = nu0 = int_list$nu0
	mean_nu1 = nu1 = int_list$nu1
	mean_logrho = log(rho)
	sigma2_logrho = int_list$sigma2_logrho
	sigma2_nu0 = int_list$sigma2_nu0
	
	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 2 + (ncol(X_nu0) + 1) + (ncol(X_logrho) + 1) + 1
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_rho = post_nu0 = post_nu1 = matrix(0, n_save, N)
	post_loglik =rep(0, n_save)
	start_time = proc.time()
	save = 1
	neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	for(k in 1:n_iter){
		
		#### Step 1 ####
		SS = SS_tv_sigma2_eps(neg_loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							  bounds=list(lower=prior_list$eps_lb, upper=prior_list$eps_ub), 
							  prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps

		SS = SS_tv_sigma2_xi(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							 bounds=list(lower=prior_list$xi_lb, upper=prior_list$xi_ub), 
							 prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
		sigma2_xi = SS$sigma2_xi
		
		# fixed time effect #
		SS = SS_tv_nu1(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		nu1 = SS$nu1	
		#### Step 2 ####
		nu0 = SS_tv_nu0_par(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mu0, C0)
		rho = SS_tv_rho_par(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0)
		
		#### step 3 ####
		betas_logrho = conj_normal(log(rho), X_logrho, sigma2_logrho, prior=list(m0=prior_list$betas0, v0=prior_list$v_betas0))
		gammas = conj_normal(nu0, X_nu0, sigma2_nu0, prior=list(m0=prior_list$gammas0, v0=prior_list$v_gammas0))
		
		#### Step 4 ####
		mean_nu0 = X_nu0%*%gammas
		res_nu0 = nu0-mean_nu0
		sigma2_nu0 = 1/rgamma(1, prior_list$a0 + (N-ncol(X_nu0))/2, prior_list$b0 + sum(res_nu0^2)/2)
		
		mean_logrho = X_logrho %*% betas_logrho
		res_logrho = log(rho) - mean_logrho
		sigma2_logrho = 1/rgamma(1, prior_list$a0 + (N-ncol(X_logrho))/2, prior_list$b0 + sum(res_logrho^2)/2)
					
		neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			post_params[save,] = c(sigma2_eps, sigma2_xi, betas_logrho, gammas, nu1[1], sigma2_logrho, sigma2_nu0)
			post_loglik[save] = 2*neg_loglik 
			post_rho[save,] = rho
			post_nu0[save,] = nu0
			post_nu1[save,] = nu1
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min

	colnames(post_params) = c("sigma2_eps", "sigma2_xi", 
							paste("betas_logrho", 0:(ncol(X_logrho)-1), sep=""),
							paste("gammas", 0:(ncol(X_nu0)-1), sep=""), 
							"nu1",
							"sigma2_logrho", "sigma2_nu0")

	summary_params = matrix(nrow = n_params, ncol=5)
	for(i in 1:n_params){
		summary_params[i,1] = mean(post_params[,i])
		summary_params[i,2] = sd(post_params[,i])
		summary_params[i,3:5] =  quantile(post_params[,i], prob = c(0.5, 0.025, 0.975))
	}
	rownames(summary_params) = colnames(post_params)
	colnames(summary_params) = c("Mean", "SD", "Med", "LCI", "UCI")

	####################################### Calculate DIC and P_d ####################################
	mean_params = colMeans(post_params)
	D_bar = mean(post_loglik)
	neg_loglik = Uni_loglik_tv(Y, time_list, mean_params[1], colMeans(post_rho), colMeans(post_nu0), colMeans(post_nu1), mean_params[2], mu0, C0)
	
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params, post_rho = post_rho, post_nu0 = post_nu0, 
				post_loglik = post_loglik))
}


######################################################    SSVM-OU ############################################################
BSSVMOU_tv2 = function(Y, X_logrho, X_nu0, X_nu1, time_list, mu0 = rep(list(c(0,0)), length(Y)), C0 = rep(list(diag(10^5, 2)),length(Y)), 
					smooth_fix = FALSE,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.1, rho = rep(5, length(Y)), 
									nu0 = rep(0, length(Y)), nu1 = rep(0, length(Y)), 
									sigma2_logrho = .1, sigma2_nu0 = .1, sigma2_nu1 = .1),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 200,
									  a0 = .1, b0 = .1,
									  betas0 = rep(0, ncol(X_logrho)), v_betas0 = rep(10^8, ncol(X_logrho)),
									  gammas0 = rep(0, ncol(X_nu0)), v_gammas0 = rep(10^8, ncol(X_nu0)),
									  etas0 = rep(0, ncol(X_nu1)), v_etas0 = rep(10^8, ncol(X_nu1))),
					n_iter = 1500, burn_in=500, thin=1, per_print=500){
	#### assign initial parameters ####
	N = length(Y)
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	mean_nu0 = nu0 = int_list$nu0
	mean_nu1 = nu1 = int_list$nu1
	mean_logrho = log(rho)
	sigma2_logrho = int_list$sigma2_logrho
	sigma2_nu0 = int_list$sigma2_nu0
	sigma2_nu1 = int_list$sigma2_nu1
	
	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 2 + (ncol(X_nu0) + 1) + (ncol(X_logrho) + 1) + (ncol(X_nu1)+1)
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_rho = post_nu0 = post_nu1 = matrix(0, n_save, N)
	post_loglik =rep(0, n_save)
	start_time = proc.time()
	save = 1
	neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
	for(k in 1:n_iter){
		
		#### Step 1 ####
		SS = SS_tv_sigma2_eps(neg_loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							  bounds=list(lower=prior_list$eps_lb, upper=prior_list$eps_ub), 
							  prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps
		if(smooth_fix == FALSE){
		SS = SS_tv_sigma2_xi(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							 bounds=list(lower=prior_list$xi_lb, upper=prior_list$xi_ub), 
							 prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
		sigma2_xi = SS$sigma2_xi
		}else{sigma2_xi = int_list$sigma2_xi}
		
		#### Step 2 ####
		nus = SS_tv_nus_par(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_nu0, sigma2_nu0, mean_nu1, sigma2_nu1, mu0, C0)
		nu0 = nus[,1]
		nu1 = nus[,2]
		rho = SS_tv_rho_par(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0)
		
		#### step 3 ####
		betas_logrho = conj_normal(log(rho), X_logrho, sigma2_logrho, prior=list(m0=prior_list$betas0, v0=prior_list$v_betas0))
		gammas = conj_normal(nu0, X_nu0, sigma2_nu0, prior=list(m0=prior_list$gammas0, v0=prior_list$v_gammas0))
		etas = conj_normal(nu1, X_nu1, sigma2_nu1, prior=list(m0=prior_list$etas0, v0=prior_list$v_etas0))

		#### Step 4 ####
		mean_nu0 = X_nu0%*%gammas
		res_nu0 = nu0-mean_nu0
		sigma2_nu0 = 1/rgamma(1, prior_list$a0 + (N-ncol(X_nu0))/2, prior_list$b0 + sum(res_nu0^2)/2)
		
		mean_nu1 = X_nu1%*%etas
		res_nu1 = nu1-mean_nu1
		sigma2_nu1 = 1/rgamma(1, prior_list$a0 + (N-ncol(X_nu1))/2, prior_list$b0 + sum(res_nu1^2)/2)
		
		mean_logrho = X_logrho %*% betas_logrho
		res_logrho = log(rho) - mean_logrho
		sigma2_logrho = 1/rgamma(1, prior_list$a0 + (N-ncol(X_logrho))/2, prior_list$b0 + sum(res_logrho^2)/2)
					
		neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			post_params[save,] = c(sigma2_eps, sigma2_xi, betas_logrho, gammas, etas, sigma2_logrho, sigma2_nu0, sigma2_nu1)
			post_loglik[save] = 2*neg_loglik 
			post_rho[save,] = rho
			post_nu0[save,] = nu0
			post_nu1[save,] = nu1
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min

	colnames(post_params) = c("sigma2_eps", "sigma2_xi", 
							paste("betas_logrho", 0:(ncol(X_logrho)-1), sep=""),
							paste("gammas", 0:(ncol(X_nu0)-1), sep=""), 
							paste("etas", 0:(ncol(X_nu1)-1), sep=""),
							"sigma2_logrho", "sigma2_nu0", "sigma2_nu1")

	summary_params = matrix(nrow = n_params, ncol=5)
	for(i in 1:n_params){
		summary_params[i,1] = mean(post_params[,i])
		summary_params[i,2] = sd(post_params[,i])
		summary_params[i,3:5] =  quantile(post_params[,i], prob = c(0.5, 0.025, 0.975))
	}
	rownames(summary_params) = colnames(post_params)
	colnames(summary_params) = c("Mean", "SD", "Med", "LCI", "UCI")

	####################################### Calculate DIC and P_d ####################################
	mean_params = colMeans(post_params)
	D_bar = mean(post_loglik)
	neg_loglik = Uni_loglik_tv(Y, time_list, mean_params[1], colMeans(post_rho), colMeans(post_nu0), colMeans(post_nu1), mean_params[2], mu0, C0)
	
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params, post_rho = post_rho, post_nu0 = post_nu0, post_nu1 = post_nu1, 
				post_loglik = post_loglik))
}



######################################################  SSVM-OU ############################################################

#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
SS_tv_nus = function(neg_loglik, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
						bounds=list(lower1=-100, upper1=100, lower2 = -100, upper2 = 100), 
						prior=list(m01=0, m02=0, v01=10^8, v02=10^8))
{
	N = length(Y_list)
	v = rexp(1) + neg_loglik
	sampled = 0; 
	lower1 = bounds$lower1; upper1 = bounds$upper1;
	lower2 = bounds$lower2; upper2 = bounds$upper2;
	while(sampled == 0){
		cand_nu0 = rtruncnorm(1, a=lower1, b=upper1, mean=prior$m01, sd=sqrt(prior$v01))
		cand_nu1 = rtruncnorm(1, a=lower2, b=upper2, mean=prior$m02, sd=sqrt(prior$v02))
		
		neg_loglik = Uni_loglik_tv(Y_list, time_list, sigma2_eps, rho, cand_nu0, cand_nu1, sigma2_xi, mu0, C0)
		
		if(v > neg_loglik){sampled = 1; nu0 = cand_nu0; nu1 = cand_nu1;}
		if( nu0 < cand_nu0){upper1 = cand_nu0}else{lower1 = cand_nu0}
		if( nu1 < cand_nu1){upper2 = cand_nu1}else{lower2 = cand_nu1}
	}
	return(list(nu0=nu0, nu1 = nu1, loglik=neg_loglik))
}


#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
SS_tv_rho = function(neg_loglik, Y_list, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
						bounds=list(lower=-7, upper=7), prior=list(m0=0, v0=10^8))
{
	N = length(Y_list)
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_rho = exp(rtruncnorm(1, a=lower, b=upper, mean=prior$m0, sd=sqrt(prior$v0)))
		neg_loglik = Uni_loglik_tv(Y_list, time_list, sigma2_eps, cand_rho, nu0, nu1, sigma2_xi, mu0, C0)
		
		if(v > neg_loglik){sampled = 1; rho = cand_rho;}
		if( log(rho) < log(cand_rho)){upper = log(cand_rho)}else{lower = log(cand_rho)}
	}
	return(list(rho=rho, loglik = neg_loglik))
}


BSSVMOU_tv_i = function(Y, time_list, mu0 = rep(list(c(0,0)), length(Y)), C0 = rep(list(diag(10^5, 2)),length(Y)), 
					smooth_fix = FALSE, lambda = 1,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.1, rho = rep(0.1, length(Y)), 
									nu0 = rep(20, length(Y)), nu1 = rep(-4, length(Y))),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 200),
					n_iter = 1500, burn_in=500, thin=1, per_print=500){
	#### assign initial parameters ####
	N = length(Y)
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	nu0 = int_list$nu0
	nu1 = int_list$nu1
	
	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 5
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_loglik =rep(0, n_save)
	start_time = proc.time()
	save = 1
	neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
									
	for(k in 1:n_iter){
		SS = SS_tv_sigma2_eps(neg_loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							  bounds=list(lower=prior_list$eps_lb, upper=prior_list$eps_ub), 
							  prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps

		SS = SS_tv_sigma2_xi(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0,
							 bounds=list(lower=prior_list$xi_lb, upper=prior_list$xi_ub), 
							 prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
		sigma2_xi = SS$sigma2_xi
		
		SS = SS_tv_nus(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		nu1 = SS$nu1
		nu0 = SS$nu0
		
		SS = SS_tv_rho(SS$loglik, Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		rho = SS$rho
		
		neg_loglik = Uni_loglik_tv(Y, time_list, sigma2_eps, rho, nu0, nu1, sigma2_xi, mu0, C0)
		c(sigma2_eps, sigma2_xi, rho, nu0, nu1)
		
		
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			post_params[save,] = c(sigma2_eps, sigma2_xi, rho, nu0, nu1)
			post_loglik[save] = 2*neg_loglik 
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min

	colnames(post_params) = c("sigma2_eps", "sigma2_xi", "rho", "nu0", "nu1")

	summary_params = matrix(nrow = n_params, ncol=5)
	for(i in 1:n_params){
		summary_params[i,1] = mean(post_params[,i])
		summary_params[i,2] = sd(post_params[,i])
		summary_params[i,3:5] =  quantile(post_params[,i], prob = c(0.5, 0.025, 0.975))
	}
	rownames(summary_params) = colnames(post_params)
	colnames(summary_params) = c("Mean", "SD", "Med", "LCI", "UCI")

	####################################### Calculate DIC and P_d ####################################
	mean_params = colMeans(post_params)
	D_bar = mean(post_loglik)
	neg_loglik = Uni_loglik_tv(Y, time_list, mean_params[1], mean_params[3], mean_params[4], mean_params[5], mean_params[2], mu0, C0)
	
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params,  
				post_loglik = post_loglik))
}

