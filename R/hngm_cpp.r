
library(truncnorm)

ngm_Y_Pred <- function(theta_old, del, rho, nu_bar, sigma2_xi, sigma2_eps){
		x1 = exp(-rho*del); x2 = exp(-2*rho*del);
		GG = matrix(0, nrow=2, ncol=3)
		GG[1,1] = 1
		GG[1,2] = (1-x1)/rho
		GG[1,3] = (del - GG[1,2])*nu_bar
		
		GG[2,2] = x1
		GG[2,3] = (1 - x1)*nu_bar

		m = GG %*% theta_old
		
		W = diag(2)
		W[1,1] = del/rho^2 + 0.5*(-3 + 4*x1 - x2)/rho^3
		W[1,2] = W[2,1] = 0.5*(1-2*x1 + x2)/rho^2
		W[2,2] = 0.5*(1-x2)/rho
		W = sigma2_xi*W
			
		theta_new = c(t(chol(W)) %*% rnorm(2) + m, 1)
		return(list(theta=theta_new, Y = rnorm(1, theta_new[1], sqrt(sigma2_eps)), mu = m, Sigma = W))
}


simul_hngm = function(n_i, sigma2_eps, sigma2_xi, rhos, nu_bars, theta0, del0 = 0.01){
	N = length(n_i)
	theta0 = cbind(theta0, 1)
	Y = theta = as.list(NULL)
	for(i in 1:N){
		theta[[i]] = matrix(0, nrow = n_i[i], ncol = 3)
		Y[[i]] = rep(0, n_i[i])
		tmp=ngm_Y_Pred(theta0[i,], del0, rhos[i], nu_bars[i], sigma2_xi, sigma2_eps)
		theta[[i]][1,] = tmp$theta
		Y[[i]][1] = tmp$Y
		for(j in 2:n_i[i]){
			tmp=ngm_Y_Pred(theta[[i]][j-1,], del_list[[i]][j], rhos[i], nu_bars[i], sigma2_xi, sigma2_eps)
			Y[[i]][j] = tmp$Y
			theta[[i]][j,] = tmp$theta
		}	
	}
	return(list(Y = Y, theta=theta))
}


dat_list = function(indata, subj_id, time_var, Y_var, del0=0.01){
	#subj_id = "subj"
	#time_var = "year"
	#Y_var = "Y"
	
	name_dat = names(indata)
	colidx_subj = which(name_dat==subj_id)
	colidx_Y = which(name_dat == Y_var)
	colidx_time = which(name_dat == time_var)
	id = unique(indata[,colidx_subj])
	N = length(id)
	find_idx_i = function(i, id_var, id) which(id_var == id[i])
	find_time_i = function(i, time_var, id_var, id) time_var[which(id_var == id[i])]
	find_Y_i = function(i, Y_var, id_var, id) Y_var[which(id_var == id[i])]
	find_del_i = function(i, time_var, id_var, id) c(del0, diff(time_var[which(id_var == id[i])]))
	idx_list = lapply(X = 1:N, find_idx_i, id_var = indata[,colidx_subj], id=id)
	time_list = lapply(X = 1:N, find_time_i, time_var = indata[,colidx_time], id_var = indata[,colidx_subj], id=id)
	del_list = lapply(X = 1:N, find_del_i, time_var = indata[,colidx_time], id_var = indata[,colidx_subj], id=id)
	
	Y_list = lapply(X = 1:N, find_Y_i, Y_var = indata[,colidx_Y], id_var = indata[,colidx_subj], id=id)
	
	return(list(idx_list = idx_list, Y_list = Y_list, time_list = time_list, del_list = del_list))
}


cpp_loglik_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				loglik = rcpp_loglik(Y_obs=Y_i, del_obs=del_list[[i]], 
									V_eps = sigma2_eps, 
									V_xi = sigma2_xi, 
									Rho=rho[i], 
									mu_int= mu00, C_int=C00)
				return(loglik)
}

loglik_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  return(sum(vapply(X=1:N, cpp_loglik_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list, FUN.VALUE = 0)))
}

cpp_theta_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				mu = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)
				return(mu)
}

cpp_thetas_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				thetas = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)[-1,]
				return(thetas)
}

SampleTheta_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  thetas = sapply(X=1:N, cpp_thetas_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
			  return(t(thetas))
}

cpp_mu_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				mu = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)[-1, 1]
				return(mu)
}

SampleMu_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  mu = sapply(X=1:N, cpp_mu_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
			  return(t(mu))
}

cpp_mu0_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				mu = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)[1, 1]
				return(mu)
}

SampleMu0_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  mu = sapply(X=1:N, cpp_mu0_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
			  return(t(mu))
}

cpp_v_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				mu = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)[-1, 2]
				return(mu)
}

SampleV_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  mu = sapply(X=1:N, cpp_v_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
			  return(t(mu))
}


cpp_v0_i = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
				if(is.list(Y)) Y_i = Y[[i]]
				if(is.matrix(Y)) Y_i = Y[i,]

				J = length(del_list[[i]])
				Z = matrix(rnorm(2*(J+1)), nrow = J+1, ncol=2)
				C00 = matrix(0, nrow = 3, ncol=3)
				C00[1:2, 1:2] = C0[[i]]; C00[3,3] = 1e-20;
				mu00 = c(mu0[[i]], nu_bar[i])
				mu = rcpp_sample_mu(Y_i, del_list[[i]], sigma2_eps, sigma2_xi, rho[i], mu00, C00, Z)[1, 2]
				return(mu)
}

SampleV0_cpp = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list){
			  if(is.list(Y)) N = length(Y)
			  if(is.matrix(Y)) N = nrow(Y)
			  mu = sapply(X=1:N, cpp_v0_i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
			  return(t(mu))
}

#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
Uni_SS_sigma2_eps = function(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi,  mu0, C0, del_list,
						bounds=list(lower=0, upper=50), prior=list(a0=.1, b0=.001))
{
	
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_sigma2_eps = rinvgamma_trunc(1, a=lower, b=upper, shape=prior$a0, rate=prior$b0)
		neg_loglik = loglik_cpp(Y, rho, nu_bar, cand_sigma2_eps, sigma2_xi, mu0, C0, del_list)
		
		if(v > neg_loglik){sampled = 1; sigma2_eps = cand_sigma2_eps;}
		if( sigma2_eps < cand_sigma2_eps){upper = cand_sigma2_eps}else{lower = cand_sigma2_eps}
	}
	return(list(sigma2_eps=sigma2_eps, loglik=neg_loglik))
}

#### shrinkage slice sampling for sigma2_xi: inverse gamma dist. was used as a prior.
Uni_SS_sigma2_xi = function(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
						bounds=list(lower=0, upper=50), prior=list(a0=0.1, b0=0.1))
{
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_sigma2_xi = rinvgamma_trunc(1, a=lower, b=upper, shape=prior$a0, rate=prior$b0)
		neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, cand_sigma2_xi, mu0, C0, del_list)
		
		if(v > neg_loglik){sampled = 1; sigma2_xi = cand_sigma2_xi;}
		if( sigma2_xi < cand_sigma2_xi){upper = cand_sigma2_xi}else{lower = cand_sigma2_xi}
	}
	return(list(sigma2_xi=sigma2_xi, loglik=neg_loglik))
}

#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
Uni_SS_rho = function(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
						bounds=list(lower=0, upper=50), prior=list(a0=.1, b0=.1))
{
   if(is.list(Y)) N = length(Y)
   if(is.matrix(Y)) N = nrow(Y)	
    current_val = rho[1]
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_val = rinvgamma_trunc(1, a=lower, b=upper, shape=prior$a0, rate=prior$b0)
		cand_rho = rep(cand_val, N)
		neg_loglik = loglik_cpp(Y, cand_rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
		
		if(v > neg_loglik){sampled = 1;}
		if( current_val < cand_val){upper = cand_val}else{lower = cand_val}
	}
	return(list(rho=cand_rho, loglik=neg_loglik))
}


#### shrinkage slice sampling for sigma2_eps: inverse gamma dist. was used as a prior.
Uni_SS_nu_bar = function(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
						bounds=list(lower=-50, upper=50), prior=list(m0=0, s0=10))
{
   if(is.list(Y)) N = length(Y)
   if(is.matrix(Y)) N = nrow(Y)	
    current_val = nu_bar[1]
	v = rexp(1) + neg_loglik
	sampled = 0; lower = bounds$lower; upper = bounds$upper;
	while(sampled == 0){
		cand_val = rtruncnorm(1, a=lower, b=upper, mean=prior$m0, sd=prior$s0)
		cand_nu_bar = rep(cand_val, N)
		neg_loglik = loglik_cpp(Y, rho, cand_nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
		
		if(v > neg_loglik){sampled = 1;}
		if( current_val < cand_val){upper = cand_val}else{lower = cand_val}
	}
	return(list(nu_bar=cand_nu_bar, loglik=neg_loglik))
}

############## sample for rho ##############
Sample_rho = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0, del_list)
{
	neg_loglik_i = cpp_loglik_i(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
	v = rexp(1) + neg_loglik_i
	cand_rho = rho

	sampled = 0; 
	lower = -7; upper =  7; 
	while(sampled == 0){
		cand_rho[i] = exp(rtruncnorm(1, a=lower, b=upper, mean=mean_logrho[i], sd=sqrt(sigma2_logrho)))

		neg_loglik_i = cpp_loglik_i(i, Y, cand_rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)

		if(v > neg_loglik_i){sampled = 1; rho[i] = cand_rho[i]; }
		if(log(rho[i]) < log(cand_rho[i])){upper = log(cand_rho[i])}else{lower = log(cand_rho[i])}
	}
	return(rho[i])
}

Sample_rho_par = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0, del_list)  
{
    if(is.list(Y)) N = length(Y)
    if(is.matrix(Y)) N = nrow(Y)
	posterior = sapply(X=1:N, Sample_rho, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0, del_list)
	return(posterior)
}

Sample_nu_bar = function(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_nu_bar, sigma2_nu_bar, mu0, C0, del_list)
{

	neg_loglik_i = cpp_loglik_i(i, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
	v = rexp(1) + neg_loglik_i 
	cand_nu_bar = nu_bar
	sampled = 0; 
	lower = -10; upper =  10; 
	while(sampled == 0){
		cand_nu_bar[i] = rtruncnorm(1, a=lower, b=upper, mean=mean_nu_bar[i], sd=sqrt(sigma2_nu_bar))
		neg_loglik_i = cpp_loglik_i(i, Y, rho, cand_nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)

		if(v > neg_loglik_i){sampled = 1; nu_bar[i] = cand_nu_bar[i];}
		if( nu_bar[i] < cand_nu_bar[i]){upper = cand_nu_bar[i]}else{lower = cand_nu_bar[i]}
	}
	return(nu_bar[i])
}


Sample_nu_bar_par = function(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_nu_bar, sigma2_nu_bar, mu0, C0, del_list)
{
    if(is.list(Y)) N = length(Y)
    if(is.matrix(Y)) N = nrow(Y)
	posterior = sapply(X=1:N, Sample_nu_bar, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_nu_bar, sigma2_nu_bar, mu0, C0, del_list)
	return(posterior)
}

######################################################    SSVM-OU ############################################################
hngm = function(Y, X_logrho, X_nu_bar, del_list, mu0 = c(0,0), C0 = diag(c(10000, 10000)), 
					smooth_fix = FALSE, lambda = 1,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.1, rho = rep(5, N), 
									nu_bar = rep(-0.5, N), sigma2_logrho = .1, sigma2_nu_bar = .1),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 200,
									  a0 = .01, b0 = .01),
					n_iter = 1500, burn_in=500, thin=1, per_print=500){
	#### assign initial parameters ####
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	nu_bar = int_list$nu_bar
	mean_logrho = log(rho)
	mean_nu_bar = nu_bar
	sigma2_logrho = int_list$sigma2_logrho
	sigma2_nu_bar = int_list$sigma2_nu_bar
	
	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 2 + (ncol(X_nu_bar) + 1) + (ncol(X_logrho) + 1) + 1
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_rho = post_nu_bar = matrix(0, n_save, N)
	post_loglik =rep(0, n_save)

	start_time = proc.time()
	save = 1
	neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
	for(k in 1:n_iter){
		#### Step 1 ####
		SS = Uni_SS_sigma2_eps(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list, 
								bounds = list(lower = prior_list$eps_lb, upper = prior_list$eps_ub), 
								 prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps
		neg_loglik = SS$loglik
		
		if(smooth_fix == FALSE){
			SS = Uni_SS_sigma2_xi(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
								  bounds = list(lower = prior_list$xi_lb, upper = prior_list$xi_ub),
								  prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
			sigma2_xi = SS$sigma2_xi
			neg_loglik = SS$loglik
		}
		#if(smooth_fix == TRUE) sigma2_xi = sigma2_eps/lambda
		if(smooth_fix == TRUE) sigma2_xi = int_list$sigma2_xi

		#### Step 2 ####
		rho = Sample_rho_par(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_logrho, sigma2_logrho, mu0, C0, del_list)	
		nu_bar = Sample_nu_bar_par(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_nu_bar, sigma2_nu_bar, mu0, C0, del_list)	

		#### Step 3 ####
		betas_logrho = conj_normal(log(rho), X_logrho, sigma2_logrho)
		betas_nu_bar = conj_normal(nu_bar, X_nu_bar, sigma2_nu_bar)
			
		#### Step 4 ####
		mean_logrho = X_logrho %*% betas_logrho
		mean_nu_bar = X_nu_bar %*% betas_nu_bar	
		
		sigma2_logrho = 1/rgamma(1, prior_list$a0 + (N-ncol(X_logrho))/2, prior_list$b0 + sum((log(rho) - mean_logrho)^2/2))
		sigma2_nu_bar = 1/rgamma(1, prior_list$a0 + (N-ncol(X_nu_bar))/2, prior_list$b0 + sum((nu_bar - mean_nu_bar)^2/2))
		
		neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
		
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			if(smooth_fix == FALSE) lambda = sigma2_eps/sigma2_xi
			post_params[save,] = c(sigma2_eps, sigma2_xi, betas_logrho, betas_nu_bar, sigma2_logrho, sigma2_nu_bar, lambda)
			post_loglik[save] = 2*neg_loglik 
			post_rho[save,] = rho
			post_nu_bar[save,] = nu_bar
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min

	colnames(post_params) = c("sigma2_eps", "sigma2_xi", 
							paste("betas_logrho", 0:(ncol(X_logrho)-1), sep=""),
							paste("betas_nu_bar", 0:(ncol(X_nu_bar)-1), sep=""),
							"sigma2_logrho", "sigma2_nu_bar", "lambda")

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
	neg_loglik = loglik_cpp(Y, colMeans(post_rho), colMeans(post_nu_bar), 
							   mean_params[1], mean_params[2], mu0, C0, del_list)
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params, post_rho = post_rho, post_nu_bar = post_nu_bar, post_loglik = post_loglik))
}


ngm = function(Y, del_list, mu0 = c(0,0), C0 = diag(c(10000, 10000)), 
					smooth_fix = FALSE, lambda = 1,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.5, rho = 5, nu_bar = -0.2),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = 1, b_xi = 1, xi_lb = 0, xi_ub = 200),
					n_iter = 1100, burn_in=100, thin=1, per_print=100){
	#### assign initial parameters ####
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	nu_bar = int_list$nu_bar

	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 4
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_loglik =rep(0, n_save)
	start_time = proc.time()
	save = 1
	neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
	for(k in 1:n_iter){
		#### Step 1 ####
		SS = Uni_SS_sigma2_eps(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list, 
								bounds = list(lower = prior_list$eps_lb, upper = prior_list$eps_ub), 
								 prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps
		neg_loglik = SS$loglik
		
		if(smooth_fix == FALSE){
			SS = Uni_SS_sigma2_xi(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
								  bounds = list(lower = prior_list$xi_lb, upper = prior_list$xi_ub),
								  prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
			sigma2_xi = SS$sigma2_xi
			neg_loglik = SS$loglik
		}
		if(smooth_fix == TRUE) sigma2_xi = sigma2_eps/lambda
		
		#### Step 2 ####
		SS = Uni_SS_rho(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)	
		rho = SS$rho
		neg_loglik = SS$loglik
		
		SS = Uni_SS_nu_bar(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)	
		nu_bar =SS$nu_bar
		neg_loglik = SS$loglik
				
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			post_params[save,] = c(sigma2_eps, sigma2_xi, rho, nu_bar)
			post_loglik[save] = 2*neg_loglik 
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min
	colnames(post_params) = c("sigma2_eps", "sigma2_xi", "rho", "nu_bar")

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
	neg_loglik = loglik_cpp(Y, mean_params[3],mean_params[4], 
							   mean_params[1], mean_params[2], mu0, C0, del_list)
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params, post_loglik = post_loglik))
}


######################################################    SSVM-Bin ############################################################

hngm_Bin = function(Y, X_nu_bar, del_list, mu0 = c(0,0), C0 = diag(c(10000, 10000)), 
					smooth_fix = FALSE, lambda = 1,
				    int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.5, rho = rep(5, N), 
									nu_bar = rep(-0.5, N), sigma2_logrho = .1, sigma2_nu_bar = .1),
				    prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
									  a_xi = 1, b_xi = 1, xi_lb = 0, xi_ub = 200,
									  a0 = .1, b0 = .1),
					n_iter = 1100, burn_in=100, thin=1, per_print=100){
	#### assign initial parameters ####
	sigma2_eps = int_list$sigma2_eps
	sigma2_xi = int_list$sigma2_xi
	rho = int_list$rho
	nu_bar = int_list$nu_bar
	mean_nu_bar = nu_bar
	sigma2_nu_bar = int_list$sigma2_nu_bar
	
	#### index for saving posterior samples ####
	idx_post = seq.int(from=(burn_in+1), to=n_iter, by=thin)
	n_save = length(idx_post)
	n_iter = tail(idx_post,1)

	#### to save posterior samples ####
	n_params = 2 + 1 + (ncol(X_nu_bar) + 1) + 1
	post_params = matrix(0, nrow = n_save, ncol=n_params)
	post_nu_bar = matrix(0, n_save, N)
	post_rho =  post_loglik =rep(0, n_save)

	start_time = proc.time()
	save = 1
	neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
	for(k in 1:n_iter){
		#### Step 1 ####
		SS = Uni_SS_sigma2_eps(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list, 
								bounds = list(lower = prior_list$eps_lb, upper = prior_list$eps_ub), 
								 prior=list(a0=prior_list$a_eps, b0=prior_list$b_eps))
		sigma2_eps = SS$sigma2_eps
		neg_loglik = SS$loglik
		
		if(smooth_fix == FALSE){
			SS = Uni_SS_sigma2_xi(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list,
								  bounds = list(lower = prior_list$xi_lb, upper = prior_list$xi_ub),
								  prior=list(a0=prior_list$a_xi, b0=prior_list$b_xi))
			sigma2_xi = SS$sigma2_xi
			neg_loglik = SS$loglik
		}
		if(smooth_fix == TRUE) sigma2_xi = sigma2_eps/lambda
		
		#### Step 2 ####
		SS = Uni_SS_rho(neg_loglik, Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
		rho = SS$rho
		
		nu_bar = Sample_nu_bar_par(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mean_nu_bar, sigma2_nu_bar, mu0, C0, del_list)	
		betas_nu_bar = conj_normal(nu_bar, X_nu_bar, sigma2_nu_bar)
		mean_nu_bar = X_nu_bar %*% betas_nu_bar	
		sigma2_nu_bar = 1/rgamma(1, prior_list$a0 + (N-ncol(X_nu_bar))/2, prior_list$b0 + sum((nu_bar - mean_nu_bar)^2/2))
		
		neg_loglik = loglik_cpp(Y, rho, nu_bar, sigma2_eps, sigma2_xi, mu0, C0, del_list)
		
		#### save sampled posterior parameters ####
		if(k == idx_post[save]){
			lambda = sigma2_eps/sigma2_xi
			post_params[save,] = c(sigma2_eps, sigma2_xi, rho[1], betas_nu_bar, sigma2_nu_bar, lambda)
			post_loglik[save] = 2*neg_loglik 
			post_rho[save] = rho[1]
			post_nu_bar[save,] = nu_bar
			save = save+1
		}
		if(0 == k %% per_print) print(paste("iteration #:", k))
	}

	end_time = proc.time()
	time_min = (end_time - start_time)[3]/60
	time_min

	colnames(post_params) = c("sigma2_eps", "sigma2_xi", "rho",
							paste("betas_nu_bar", 0:(ncol(X_nu_bar)-1), sep=""),
							"sigma2_nu_bar", "lambda")

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
	neg_loglik = loglik_cpp(Y, rep(mean(post_rho), N), colMeans(post_nu_bar), 
							   mean_params[1], mean_params[2], mu0, C0, del_list)
	D_hat = 2*neg_loglik
	DIC = 2*D_bar - D_hat
	P_d = D_bar - D_hat
	
	return(list(summary_params = summary_params, DIC=DIC, P_d = P_d, 
				post_params = post_params, post_rho = post_rho, post_nu_bar = post_nu_bar, post_loglik = post_loglik))
}




## indicator whether poster mean/velocity is calculated with covariates adjusted
post_dynamics = function(hngm_object, X_logrho, X_nu_bar, 
						 covariate_adjust = FALSE, 
						 Y_obs, time_obs, time_aug, id = c(1:length(Y_obs)),
						 visual = FALSE, dir_figure)
{
	N = nrow(X_nu_bar)
	#### procedure to get posterior samples with adjusting X at their means ####
	X_logrho_mean =  matrix(colMeans(X_logrho), nrow = N, ncol=ncol(X_logrho), byrow = TRUE)
	X_nu_bar_mean = matrix(colMeans(X_nu_bar), nrow = N, ncol=ncol(X_nu_bar), byrow = TRUE)
	p_logrho = ncol(X_logrho)
	p_nu_bar = ncol(X_nu_bar)
	
	#### index for betas of logrho ####
	idx_logrho = 3:(2+p_logrho)
	#### index for gammas of nu_bar ####
	idx_nu_bar = (3+p_logrho):(2+ p_logrho + p_nu_bar)
	
	post_rho = hngm_object$post_rho
	post_mean_logrho = hngm_object$post_params[,idx_logrho] %*% t(X_logrho)
	b0_logrho = log(post_rho) - post_mean_logrho
	
	post_nu_bar = hngm_object$post_nu_bar
	post_mean_nu_bar = hngm_object$post_params[,idx_nu_bar] %*% t(X_nu_bar)
	b0_nu_bar = post_nu_bar - post_mean_nu_bar

	post_rho_adj = exp(hngm_object$post_params[,idx_logrho] %*% t(X_logrho_mean) + b0_logrho)
	post_nu_bar_adj = hngm_object$post_params[,idx_nu_bar] %*% t(X_nu_bar_mean) + b0_nu_bar

	post_sigma2_eps = hngm_object$post_params[,1]
	post_sigma2_xi = hngm_object$post_params[,2]
	n_post = nrow(hngm_object$post_params)

	###################################### get smoothed mean profile #######################################
	#### set time sequence to predict outcomes at those points #####
	J = length(time_aug)
	#### set missing indicators at multiple visit points ####
	miss_mat = matrix(1, nrow = N, ncol=J)
	for(i in 1:N){
		rounded_time = round(time_obs[[i]], 1)
		idx_obs = which(round(time_aug,1) %in% rounded_time)
		miss_mat[i, idx_obs] = 0
	}

	#### Set missing values of BMI with initial value (=15) for data augmentation ####
	Y_aug = matrix(0, nrow = N, ncol=J)
	for(i in 1:N){
		idx_obs = which(miss_mat[i,]==0)
		idx_miss = which(miss_mat[i,] == 1)
		Y_aug[i, idx_obs] = Y_obs[[i]]
		Y_aug[i, idx_miss] = mean(Y_obs[[i]])
	}

	#### Perform initial augmentation for missing values ####
	del_aug = rep(list(c(0.01, diff(time_aug))), N)
	N_augment = 1000
	for(k in 1:N_augment){
		if(covariate_adjust==TRUE) mu = SampleMu_cpp(Y_aug, post_rho_adj[1,], post_nu_bar_adj[1,], post_sigma2_eps[1], post_sigma2_xi[1], mu0, C0, del_aug)
		if(covariate_adjust==FALSE) mu = SampleMu_cpp(Y_aug, post_rho[1,], post_nu_bar[1,], post_sigma2_eps[1], post_sigma2_xi[1], mu0, C0, del_aug)
		Y_sample = mu + rnorm(N*J, 0, sd=sqrt(post_sigma2_eps[1]))
		Y_aug[miss_mat==1] = Y_sample[miss_mat==1]
	}

	#### Posterior sampling for mean and velocity across subjects ####
	idx_U = 1:J; idx_V = idx_U+J
	post_U = post_V = array(0, dim=c(n_post, J, N))
	for(i in 1:n_post){	
		if(covariate_adjust==TRUE) thetas = SampleTheta_cpp(Y_aug, post_rho_adj[i,], post_nu_bar_adj[i,], post_sigma2_eps[i], post_sigma2_xi[i], mu0, C0, del_aug)
		if(covariate_adjust==FALSE)thetas = SampleTheta_cpp(Y_aug, post_rho[i,], post_nu_bar[i,], post_sigma2_eps[i], post_sigma2_xi[i], mu0, C0, del_aug)
		
		post_U[i,,] = t(thetas[, idx_U])
		post_V[i,,] = t(thetas[, idx_V])
		Y_sample = thetas[, idx_U] + rnorm(N*J, 0, sd=sqrt(post_sigma2_eps[i]))
		Y_aug[miss_mat==1] = Y_sample[miss_mat==1]	
	}

	#### summarize estimates for individual mean and velocity #### 
	summary_U_i = summary_V_i = array(dim=c(length(time_aug),3, N))
	for(s in 1:N)
		for(j in 1:J){
			summary_U_i[j,,s] = quantile(post_U[,j,s], prob = c(0.5, 0.025, 0.975))
			summary_V_i[j,,s] = quantile(post_V[,j,s], prob = c(0.5, 0.025, 0.975))
		}

	if(visual==TRUE){
		######################### individual trajectory/velocity with infant BMI peak #######################
		pdf(dir_figure,width=7,height=160)
		par(mfrow = c(ceiling((N+1)/3),3),mar= c(4, 3.5, 4, 3.2) + 0.1)
		for(s in 1:N){
			range_Y = round(range(unlist(Y_obs)))
			plot(0, ylab = "", xlab = "Year", type="n", xlim=range(time_aug), ylim=range_Y, yaxt="n")
			points(time_obs[[s]], Y_obs[[s]], col="red")
			u_at = round(seq(range_Y[1], range_Y[2], by=diff(range_Y)/3),1)
			axis(2, at=u_at,labels=u_at, las=2, cex.axis=1, tck=-.02, hadj=.5)
			mtext(expression(hat(U)(t)), side=2, cex.lab=.5, las=1, adj=2, cex=.7)
			lines(time_aug, summary_U_i[,1,s])
			lines(time_aug, summary_U_i[,2,s], lty=2)
			lines(time_aug, summary_U_i[,3,s], lty=2)
			idx_v0_neg = which((summary_V_i[,1,s]<0) ==1)
			title(paste("Individual ID: ", id[s], sep=""))
			if(length(idx_v0_neg) > 0){
				idx_peak = idx_v0_neg[1]
				points(time_aug[idx_peak], summary_U_i[idx_peak,1,s], pch=8)
			}
			lgd = c("Observed BMI", "Infant peak")
			if(s==1)legend("topright", lgd, pch=c(1,8), col=c("red", "black"), bty="n")
		}
		dev.off()
	}
	return(list(summary_U_i = summary_U_i, summary_V_i = summary_V_i))
}




