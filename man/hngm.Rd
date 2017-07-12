\name{SSVM_OU}
\alias{SSVM_OU}
\docType{package}
\title{
Bayesian Semi-parametric Stochastic Model with Ornstein-Uhlenbeck process prior (B-SSVM-OU)
}
\description{
To fit B-SSVM-OU with longitudial data to examine the association between growth acceleration and an exposure.
}
\usage{
SSVM_OU(Y_list, X_logrho, X_nu_bar, del_list, mu0, C0, 
		smooth_fix = FALSE, lambda = 1,
		int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.1, rho = rep(5, N), 
						nu_bar = rep(-0.5, N), sigma2_logrho = .1, sigma2_nu_bar = .1),
		prior_list = list(a_eps = .001, b_eps = .001, eps_lb = 0, eps_ub = 10,
						  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 200, a0 = .01, b0 = .01),
		n_iter = 1500, burn_in=500, thin=1, per_print=500)
}
\examples{
\dontrun{

data(sparse_dat) ## visit for each subject J_i = 8 ##
## same for data(moderate_dat): visit for each subject J_i = 30 ##
head(sparse_dat)

dat = dat_list(indata=sparse_dat, subj_id="subj", time_var="year", Y_var="Y")	
del_list = dat$del_list
Y_list = dat$Y_list
N = length(Y_list)
X_dat = unique(data.frame(subj=sparse_dat$subj, X=sparse_dat$X))

X_logrho = model.matrix(~ X, data=X_dat)
X_nu_bar = model.matrix(~ 1, data=X_dat)

####################################    MCMC simulation   #####################################
n_iter = 2500; burn_in = 500; thin=1; per_print=500
## prior for initial distribution of trajectory and velocity (should be set as a list format) 
mu0 = C0 = as.list(NULL)
for(i in 1:N){
	mu0[[i]] = c(13.1, 31) # c(m01, m02), where m01 and m02 are prior means, respectively, for trajectory and velocity at the initial time.
	C0[[i]] = diag(c(10, 20)) # c(c01, c02), where m01 and m02 are prior variances, respectively, for trajectory and velocity at the initial time.
}

## initial values for MCMC iterations.
int_list = list(sigma2_eps = 0.5, sigma2_xi = 0.5, rho = rep(5, N), nu_bar = rep(-0.5, N),
		sigma2_logrho = .1, sigma2_nu_bar = .1)

## values for hyper-parameters.
## sigma2_eps ~ IG(a_eps, b_eps)I(eps_lb, eps_ub)
## sigma2_xi ~ IG(a_xi, b_xi)I(xi_lb, xi_ub)
## sigma2_logrho, sigma2_gamma ~ IG(a0, b0)I(0, inf).
prior_list = list(a_eps = .01, b_eps = .01, eps_lb = 0, eps_ub = 10000,
				  a_xi = .1, b_xi = .1, xi_lb = 0, xi_ub = 10000,
				  a0 = .01, b0 = .01)

SSVM = SSVM_OU(Y_list, X_logrho, X_nu_bar, del_list, mu0=mu0, C0=C0, 
			   prior_list = prior_list, int_list = int_list, 
			   n_iter = n_iter, burn_in=burn_in, thin=thin, per_print=per_print)
}
}
