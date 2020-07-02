holdout_data = function(X, holdout_portion=0.25){
	num_datapoints = dim(X)[1]
	data_dim = dim(X)[2]

    n_holdout = as.integer(holdout_portion * num_datapoints * data_dim)

    holdout_row = sample(1:num_datapoints, n_holdout, replace=TRUE)
    holdout_col = sample(1:data_dim, n_holdout, replace=TRUE)
    holdout_mask = as.matrix(sparseMatrix(i=holdout_row, j=holdout_col, 
    	x=rep(1, n_holdout), dims=c(num_datapoints, data_dim)))

    # prevent an entry from being selected twice
    holdout_mask = pmin(holdout_mask, matrix(1, num_datapoints, data_dim))

    holdout_subjects = unique(holdout_row)

    # we hold out entries by setting the heldout entries to be zero
    # moreover, we will ignore the likelihood of heldout entries in
    # fitting factor models. so setting heldout entries as zero shall
    # not bias model fits.

    x_train = (1-holdout_mask) * X
    x_vad = holdout_mask * X

    returnvals = list("x_train"=x_train, "x_vad"=x_vad, 
    	"holdout_row"=holdout_row, "holdout_col"=holdout_col, 
    	"holdout_mask"=holdout_mask)

    return(returnvals)
}



compute_beta_acc <- function(beta_fit, betas){
	# beta_fit is a matrix of beta samples \hat{\beta}

	# betas is true beta

	# consider the estimate \hat{\beta} as a randomly drawn
	# approximate posterior sample of the beta parameter

	# bias =  E[\hat{\beta}] - beta
	# var = Var(\hat{\beta}) 
	# mse = E[(\hat{\beta} - beta)^2]

	# these definitions follow the bias / variance / mse of posterior
	# samples (Korattikara et al., 2014; Chen et al., 2015).

	# Korattikara, A., Chen, Y., & Welling, M. (2014). Austerity in
	# MCMC land: Cutting the Metropolis-Hastings budget. In
	# International Conference on Machine Learning (pp. 181-189).

	# Chen, C., Ding, N., & Carin, L. (2015). On the convergence of
	# stochastic gradient MCMC algorithms with high-order integrators.
	# In Advances in Neural Information Processing Systems (pp.
	# 2278-2286).

	bias2 = sum((apply(beta_fit, 2, mean) - betas)^2)
	var = sum(apply(beta_fit, 2, var))
	mse = sum(apply((beta_fit - matrix(betas, nrow=dim(beta_fit)[1], 
		ncol=dim(beta_fit)[2], byrow=TRUE))^2, 2, mean))
	returnvals = list("bias2"=bias2, "var"=var, "mse"=mse)

	return(returnvals)
}



save_res <- function(res, config){
	setname = deparse(substitute(res))
	filename = paste(c(config, setname, ".csv"), collapse="_")
	write.table(res, file = filename, sep = ",", 
		col.names = c("bias2", "var", "mse"),
	    qmethod = "double")
}

save_simdat <- function(var, savedir, randseed, simset){
	setname = deparse(substitute(var))
	filename = paste(c(savedir, paste(c(randseed, simset, 
		setname, ".csv"), collapse="_")), collapse="/")
	write.table(var, file = filename, sep = ",", 
		qmethod = "double")
}

load_simdat <- function(varname, simdatdir, randseed, simset){
	filename = paste(c(simdatdir, paste(c(randseed, simset, 
		varname, ".csv"), collapse="_")), collapse="/")
	varvals = read.table(filename, sep=",", header=TRUE)
	return(varvals)
}

sim_conf_data <- function(A, C, dir, family=gaussian(), 
	sd=1.0, confscale=2.0, seed=0){

	# family: gaussian() or binomial() 
	# simcount: counts of this is xxth simulated dataset

	N = dim(A)[1]
	D = dim(A)[2]
	Kc = dim(C)[2]

	intercept = rnorm(n=1)
	betas = rnorm(n=D)
	gammas = rnorm(n=Kc) * confscale

	simmean = intercept + A %*% betas + C %*% gammas

	if (identical(family, gaussian(), ignore.environment=TRUE)){
		simY = rnorm(n=N, mean=simmean, sd=sd)
	} else if (identical(family, binomial(), ignore.environment=TRUE)){
		prob = 1/(1+exp(-simmean)) 
		simY = rbinom(n=N, size=1, prob=prob)   
	}

	returnvals = list("intercept"=intercept, "betas"=betas, "gammas"=gammas, "simmean"=simmean, "simY"=simY)

	# stan_rdump(c("intercept", "betas", "gammas", "simmean", "simY"), paste(dir, paste(paste("simdata", seed, sep=""), ".rdump", sep=""), sep="/"))

	#can use read_rdump to load this file

	return(returnvals)
}



