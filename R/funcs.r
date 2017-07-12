
#### See R Programs for Truncated Distributions by Saralees Nadarajah and Samuel Kotz ####
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}

rgamma_trunc = function(n, a=-Inf, b=Inf, shape, rate){
			return(rtrunc(n, "gamma", a, b, shape=shape, rate=rate))
}

rinvgamma_trunc = function(n, a=-Inf, b=Inf, shape, rate){
			if(a==0) a=1e-30
			return(1/rtrunc(n, "gamma", a=1/b, b=1/a, shape=shape, rate=rate))
}


r_trunc_gamma = function(n_sample, a, b, alpha, beta){
                     
            ## to sample truncated gamma in a and b ##
            h = function(u, a, b, alpha, beta){
            		m = min(b, -log(u)/beta)
            		val = m^alpha - a^alpha
            		return(val)
            }
            
            iCDF = function(a, b, alpha){
            		v = runif(1)
            		x = b*(v + (a/b)^alpha*(1-v))^(1/alpha)
            		return(x)
            }

     if(b>=100){ 

            X = rep(0, n_sample)
            for(j in 1:n_sample){
            
                y = runif(1, 0, exp(-beta*a))
                L_b = min(b, -log(y)/beta)
                x = iCDF(a, L_b, alpha)
                n_iter = 100
                for(i in 1:n_iter){
                       y = runif(1, 0, exp(-x*beta))
                       L_b = min(b, -log(y)/beta)
                       x = iCDF(a, L_b, alpha)
                }
                X[j] = x
            }
     }

     if(b<100){             
            L = 0; U = exp(-a*beta)
            X = rep(NA, n_sample)
            
            for(s in 1:n_sample){
            
                interval = c(L, U)
                sampled_idx = 1
                sampled = 0
                while(sampled == 0){
                
                	sampled_int = c(interval[sampled_idx], interval[sampled_idx+1])
                	u = runif(1, sampled_int[1], sampled_int[2])
                	h_L = h(sampled_int[1], a, b, alpha, beta)
                	h_u = h(u, a, b, alpha, beta)
                
                	w = runif(1)
                
                	if(w < h_u/h_L){sampled = 1}
                	if(w > h_u/h_L){
                    		interval = sort(c(interval, u))
                    		n_int = length(interval) - 1
                    		gap_int = rep(0, n_int)
                    		for(i in 1:n_int) gap_int[i] = interval[i+1] - interval[i]
                    		h_int = rep(0, n_int)
                    		for(i in 1:n_int) h_int[i] = h(interval[i], a,b, alpha, beta)
                    		prob_int = (h_int*gap_int)/sum(h_int*gap_int)
                    		sampled_idx = sample(1:n_int, size = 1, prob = prob_int)
                	}
                }
                
                L_a = a
                L_b = min(b, -log(u)/beta)
                x = iCDF(L_a, L_b, alpha)
                
                X[s] = x
            }
	    
      }
      return(X)
}


r_trunc_igamma = function(n_sample, a, b, alpha, beta){
			if(a == 0) a = 1e-30
			X = 1/r_trunc_gamma(n_sample, a = 1/b, b = 1/a, alpha=alpha, beta=beta)
		        return(X)
		 }


#### conjugate form of normal distribution to regress nu_bar or logrho_bar to their covariates.
conj_normal = function(Y, X, sigma_sq, prior = list(m0=rep(0, ncol(X)), v0=rep(10^8, ncol(X)))){
				p = ncol(X)
				Sigma = solve(t(X) %*% X/sigma_sq + diag(1/prior$v0,p))
				mu = Sigma %*% ( (t(X) %*% Y)/sigma_sq + prior$m0/prior$v0)
				mu_til =  t(chol(Sigma)) %*% rnorm(p) + mu
				return(mu_til)
}



