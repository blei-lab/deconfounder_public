rm(list = ls())

rawdir = "./smoking_R"

setwd(paste(c(rawdir, "src"), collapse="/"))

dir.create(paste(c(rawdir, "res/factorfit"), collapse="/"))


###########################################################
# parse arguments
###########################################################


args <- commandArgs(trailingOnly = TRUE)

factor_model = args[1]
factor_model_K = as.integer(args[2])
estctrl = args[3]
trial = as.integer(args[4])
simset = args[5]
dataseed = as.integer(args[6])
factorseed = as.integer(args[7])
randseed = as.integer(args[8])
algorithm = args[9]
priorsp = args[10]
outcome_model = args[11]

# factor_model = "quadratic"
# factor_model_K = 1
# estctrl = "Actrl"
# trial = 49
# simset = "indep"
# dataseed = 20133224
# factorseed = 20144133
# randseed = 20200606
# algorithm = "meanfield"
# priorsp = "normal"
# outcome_model = "linear"


# number of post sample draws 
maxsmps = 998
npostsmps = 5 
# can increase the number of posterior samples especially if needs to
# compute credible intervals


savedir = paste(c(rawdir, "res", "causalest", simset,
	paste(c(factor_model, factor_model_K), collapse="_")),
	collapse="/")

dir.create(paste(c(rawdir, "res", "causalest"),
	collapse="/"))

dir.create(paste(c(rawdir, "res", "causalest", simset),
	collapse="/"))

dir.create(savedir)

params <- data.frame(
	"factor_model" = c(factor_model), 
	"factor_model_K" = c(factor_model_K), 
	"estctrl" = c(estctrl),
	"trial" = c(trial),
	"simset" = c(simset),
	"dataseed" = c(dataseed),
	"factorseed" = c(factorseed),
	"randseed" = c(randseed),
	"algorithm" = c(algorithm), 
	"priorsp" = c(priorsp), 
	"outcome_model" = c(outcome_model))

params_filenames = paste(c(randseed, dataseed, factorseed, estctrl, "causalest_params.csv"), collapse="_")
write.csv(params, paste(c(savedir, params_filenames), collapse="/"))


###########################################################
# set up 
###########################################################

library("rstan")
library("rstanarm")
library("Matrix")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("utils.R")

randseed = as.integer(as.numeric(Sys.time()))
set.seed(randseed)
print(randseed)


###########################################################
# set up simulation results placeholders
###########################################################

default_val = NA

if(simset=="indep"){
	simdatdir = paste(c(rawdir, "dat/simdat/indep"), collapse="/")
} else if(simset=="dep"){
	# for causes with dependence	
	simdatdir = paste(c(rawdir, "dat/simdat/dep"), collapse="/")
}

A = as.matrix(load_simdat("causes", simdatdir, dataseed, simset))
C = as.matrix(load_simdat("confounders", simdatdir, dataseed, simset))
cov = as.matrix(load_simdat("covariates", simdatdir, dataseed, simset))
simYs = load_simdat("simYs", simdatdir, dataseed, simset)
betass = load_simdat("betass", simdatdir, dataseed, simset)

N = dim(A)[1]
D = dim(A)[2]
Kc = dim(C)[2]

if (priorsp=="normal"){
	prior = normal()
} else if (priorsp=="horseshoe"){
	prior = hs_plus()
}

if(outcome_model=="linear"){
	family = gaussian()
} else if (outcome_model=="logistic"){
	family = binomial()
}

simY = c(t(simYs[trial,]))
betas = c(t(betass[trial,]))


fitfactordir = paste(c(rawdir, "res", "factorfit", simset,
	paste(c(factor_model, factor_model_K), collapse="_")),
	collapse="/")

# ###########################################################
# # no confounding adjustment
# ###########################################################


if (estctrl=="nctrl"){
	nctrlfit = stan_glm.fit(cbind(rep(1,N), A), 
		simY, algorithm=algorithm, prior=prior, family=family, 
		QR=TRUE,prior_intercept = normal(),
		output_samples=npostsmps)

	la = extract(nctrlfit)
	beta_fit = la$beta[1:npostsmps,2:(D+1)]
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	nctrl = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)
	
	# for sanity check
	# fit1 = glm(simY~A, family=family)
	# summary(fit1)

	# freqnctrl_bias2 = sum((summary(fit1)$coefficients[2:(2+D-1),1] - betas)^2)
	# freqnctrl_var = sum((summary(fit1)$coefficients[2:(2+D-1),2])^2)
	# freqnctrl_mse = freqnctrl_bias2 + freqnctrl_var
	# freqnctrl = c(freqnctrl_bias2, freqnctrl_var, freqnctrl_mse)	

	# freqnctrl are not used for reporting, but are used only to cross
	# check if Bayesian regression has been properly fitted

	# nctrlres = data.frame(rbind(nctrl, freqnctrl))
	nctrlres = data.frame(rbind(nctrl))
	colnames(nctrlres) = c("bias2", "var", "mse", "cov1", "cov2")
	print(nctrlres)

	res = nctrlres
}

	



# ###########################################################
# # adjust for true confounder
# ###########################################################

if (estctrl=="oracle"){
	oraclefit = stan_glm.fit(cbind(rep(1,N), A, C), 
		simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
		prior_intercept = normal(),
		output_samples=npostsmps)

	la = extract(oraclefit)
	beta_fit = la$beta[,2:(D+1)]	
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	oracle = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)

	# for sanity check
	# fit2 = glm(simY~A+C, family=family)
	# summary(fit2)
	# freqoracle_bias2 = sum((summary(fit2)$coefficients[2:(2+D-1),1] - betas)^2)
	# freqoracle_var = sum((summary(fit2)$coefficients[2:(2+D-1),2])^2)
	# freqoracle_mse = freqoracle_bias2 + freqoracle_var
	# freqoracle = c(freqoracle_bias2, freqoracle_var, freqoracle_mse)

	
	# oracleres = data.frame(rbind(oracle, freqoracle))

	oracleres = data.frame(rbind(oracle))
	colnames(oracleres) = c("bias2", "var", "mse", "cov1", "cov2")
	print(oracleres)

	res = oracleres

}





# ###########################################################
# # adjust for substitute confounder
# ###########################################################

if (estctrl=="Zctrl"){
	beta_fit = matrix(default_val, npostsmps, D)
	smps = sample.int(maxsmps, npostsmps)
	for (j in (1:npostsmps)){
		filenameZ = paste(c(dataseed, factorseed, "Zhat", smps[j], ".csv"), collapse="_")
		Zhat = as.matrix(read.table(paste(c(fitfactordir, filenameZ), collapse="/"), 
			sep=",", header=TRUE))

		tryCatch({Zctrlfit = stan_glm.fit(cbind(rep(1,N), A, Zhat), 
					simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
					prior_intercept = normal(),
					output_samples=100)},
		error=function(e){
			cat("ERROR :",conditionMessage(e), "\n");
			Zctrlfit = stan_glm.fit(cbind(rep(1,N), A, Zhat), 
			simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
			prior_intercept = normal(),
			output_samples=100);
		})


		la = extract(Zctrlfit)
		beta_fit[j,] = la$beta[sample.int(100,1),2:(D+1)]
	}
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	Zctrl = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)

	# freq_beta_fit = matrix(default_val, npostsmps, D)
	# freq_beta_var = matrix(default_val, npostsmps, D)
	# for (j in 1:npostsmps){
	# 	filenameZ = paste(c(dataseed, factorseed, "Zhat", smps[j], ".csv"), collapse="_")
	# 	Zhat = as.matrix(read.table(paste(c(fitfactordir, filenameZ), collapse="/"), 
	# 		sep=",", header=TRUE))
	# 	fit3 = glm(simY~A+Zhat, family=family)
	# 	summary(fit3)
	# 	freq_beta_fit[j,] = summary(fit3)$coefficients[2:(2+D-1),1]
	# 	freq_beta_var[j,] = (summary(fit3)$coefficients[2:(2+D-1),2])^2
	# }
	# freqZctrl_bias2 = sum((apply(freq_beta_fit, 2, mean) - betas)^2)
	# freqZctrl_var = sum(apply(freq_beta_fit, 2, var)) + sum(apply(freq_beta_var, 2, mean))
	# freqZctrl_mse = freqZctrl_bias2 + freqZctrl_var
	# freqZctrl = c(freqZctrl_bias2, freqZctrl_var, freqZctrl_mse)

	# Zctrlres = data.frame(rbind(Zctrl, freqZctrl))

	Zctrlres = data.frame(rbind(Zctrl))
	colnames(Zctrlres) = c("bias2", "var", "mse", "cov1", "cov2")

	print(Zctrlres)

	res = Zctrlres

}





###########################################################
# adjust for reconstructed causes
###########################################################
# better use prior = hs_plus(slab_scale=0.5) for Ahat
# prior = hs_plus()

if (estctrl=="Actrl"){
	beta_fit = matrix(default_val, npostsmps, D)
	smps = sample.int(maxsmps, npostsmps)
	for (j in (1:npostsmps)){
		filenameA = paste(c(dataseed, factorseed, "Ahat", smps[j], ".csv"), collapse="_")
		Ahat = as.matrix(read.table(paste(c(fitfactordir, filenameA), collapse="/"), 
			sep=",", header=TRUE))
		tryCatch({Actrlfit = stan_glm.fit(cbind(rep(1,N), A, Ahat), 
				simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
					prior_intercept = normal(),
					output_samples=100)},
		error=function(e){
			cat("ERROR :",conditionMessage(e), "\n");
			Actrlfit = stan_glm.fit(cbind(rep(1,N), A, Ahat), 
			simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
			prior_intercept = normal(),
			output_samples=100);
		})
		la = extract(Actrlfit)
		beta_fit[j,] = la$beta[sample.int(100,1),2:(D+1)]
	}
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	Actrl = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)


	# freq_beta_fit = matrix(default_val, npostsmps, D)
	# freq_beta_var = matrix(default_val, npostsmps, D)
	# for (j in 1:npostsmps){
	# 	filenameA = paste(c(dataseed, factorseed, "Ahat", smps[j], ".csv"), collapse="_")
	# 	Ahat = as.matrix(read.table(paste(c(fitfactordir, filenameA), collapse="/"), 
	# 		sep=",", header=TRUE))
	# 	Ares = A - Ahat
	# 	fit4 = glm(simY~Ares, family=family)
	# 	summary(fit4)
	# 	freq_beta_fit[j,] = summary(fit4)$coefficients[2:(2+D-1),1]
	# 	freq_beta_var[j,] = (summary(fit4)$coefficients[2:(2+D-1),2])^2
	# }
	# freqActrl_bias2 = sum((apply(freq_beta_fit, 2, mean) - betas)^2)
	# freqActrl_var = sum(apply(freq_beta_fit, 2, var)) + sum(apply(freq_beta_var, 2, mean))
	# freqActrl_mse = freqActrl_bias2 + freqActrl_var
	# freqActrl = c(freqActrl_bias2, freqActrl_var, freqActrl_mse)

	# Actrlres = data.frame(rbind(Actrl, freqActrl))

	Actrlres = data.frame(rbind(Actrl))
	colnames(Actrlres) = c("bias2", "var", "mse", "cov1", "cov2")

	print(Actrlres)

	res = Actrlres
}




###########################################################
# adjust for substitute confounder and all covariates
###########################################################

if (estctrl=="Zcovctrl"){

	beta_fit = matrix(default_val, npostsmps, D)
	smps = sample.int(maxsmps, npostsmps)
	for (j in 1:npostsmps){
		filenameZ = paste(c(dataseed, factorseed, "Zhat", smps[j], ".csv"), collapse="_")
		Zhat = as.matrix(read.table(paste(c(fitfactordir, filenameZ), collapse="/"), 
			sep=",", header=TRUE))

		tryCatch({Zcovctrlfit = stan_glm.fit(cbind(rep(1,N), A, Zhat, cov), 
					simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family,
					prior_intercept = normal(),
					output_samples=100)},
		error=function(e){
			cat("ERROR :",conditionMessage(e), "\n");
			Zcovctrlfit = stan_glm.fit(cbind(rep(1,N), A, Zhat, cov), 
			simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family,
			prior_intercept = normal(),
			output_samples=100);
		})
		la = extract(Zcovctrlfit)
		beta_fit[j,] = la$beta[sample.int(100,1),2:(D+1)]
	}
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	Zcovctrl = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)

	# freq_beta_fit = matrix(default_val, npostsmps, D)
	# freq_beta_var = matrix(default_val, npostsmps, D)
	# for (j in 1:npostsmps){
	# 	filenameZ = paste(c(dataseed, factorseed, "Zhat", smps[j], ".csv"), collapse="_")
	# 	Zhat = as.matrix(read.table(paste(c(fitfactordir, filenameZ), collapse="/"), 
	# 		sep=",", header=TRUE))
	# 	fit3 = glm(simY~A+Zhat+cov, family=family)
	# 	summary(fit3)
	# 	freq_beta_fit[j,] = summary(fit3)$coefficients[2:(2+D-1),1]
	# 	freq_beta_var[j,] = (summary(fit3)$coefficients[2:(2+D-1),2])^2
	# }
	# freqZcovctrl_bias2 = sum((apply(freq_beta_fit, 2, mean) - betas)^2)
	# freqZcovctrl_var = sum(apply(freq_beta_fit, 2, var)) + sum(apply(freq_beta_var, 2, mean))
	# freqZcovctrl_mse = freqZcovctrl_bias2 + freqZcovctrl_var
	# freqZcovctrl = c(freqZcovctrl_bias2, freqZcovctrl_var, freqZcovctrl_mse)

	# Zcovctrlres = data.frame(rbind(Zcovctrl, freqZcovctrl))

	Zcovctrlres = data.frame(rbind(Zcovctrl))
	colnames(Zcovctrlres) = c("bias2", "var", "mse", "cov1", "cov2")
	print(Zcovctrlres)

	res = Zcovctrlres
}





###########################################################
# adjust for reconstructed causes and all covariates
###########################################################

# prior = hs_plus(slab_scale=0.5)
# prior = hs_plus()


if (estctrl=="Acovctrl"){
	beta_fit = matrix(default_val, npostsmps, D)
	smps = sample.int(maxsmps, npostsmps)
	for (j in (1:npostsmps)){
		filenameA = paste(c(dataseed, factorseed, "Ahat", smps[j], ".csv"), collapse="_")
		Ahat = as.matrix(read.table(paste(c(fitfactordir, filenameA), collapse="/"), 
			sep=",", header=TRUE))
		tryCatch({Acovctrlfit = stan_glm.fit(cbind(rep(1,N), A-Ahat, cov), 
				simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
					prior_intercept = normal(),
					output_samples=100)},
		error=function(e){
			cat("ERROR :",conditionMessage(e), "\n");
			Acovctrlfit = stan_glm.fit(cbind(rep(1,N), A-Ahat, cov), 
			simY, algorithm=algorithm, QR=TRUE, prior=prior, family=family, 
			prior_intercept = normal(),
			output_samples=100);
		})
		la = extract(Acovctrlfit)
		beta_fit[j,] = la$beta[sample.int(100,1),2:(D+1)]
	}
	beta_fit_acc = compute_beta_acc(beta_fit, betas)
	credible_int_1 = quantile(beta_fit[,1], probs=c(0.025, 0.975))
	credible_int_2 = quantile(beta_fit[,2], probs=c(0.025, 0.975))
	cov_1 = (betas[1] < credible_int_1[2]) & (betas[1] > credible_int_1[1])
	cov_2 = (betas[2] < credible_int_2[2]) & (betas[2] > credible_int_2[1])
	Acovctrl = c(beta_fit_acc$bias2, beta_fit_acc$var, beta_fit_acc$mse, cov_1, cov_2)


	# freq_beta_fit = matrix(default_val, npostsmps, D)
	# freq_beta_var = matrix(default_val, npostsmps, D)
	# for (j in 1:npostsmps){
	# 	filenameA = paste(c(dataseed, factorseed, "Ahat", smps[j], ".csv"), collapse="_")
	# 	Ahat = as.matrix(read.table(paste(c(fitfactordir, filenameA), collapse="/"), 
	# 		sep=",", header=TRUE))
	# 	fit4 = glm(simY~A-Ahat+cov, family=family)
	# 	summary(fit4)
	# 	freq_beta_fit[j,] = summary(fit4)$coefficients[2:(2+D-1),1]
	# 	freq_beta_var[j,] = (summary(fit4)$coefficients[2:(2+D-1),2])^2
	# }
	# freqAcovctrl_bias2 = sum((apply(freq_beta_fit, 2, mean) - betas)^2)
	# freqAcovctrl_var = sum(apply(freq_beta_fit, 2, var)) + sum(apply(freq_beta_var, 2, mean))
	# freqAcovctrl_mse = freqAcovctrl_bias2 + freqAcovctrl_var
	# freqAcovctrl = c(freqAcovctrl_bias2, freqAcovctrl_var, freqAcovctrl_mse)

	# Acovctrlres = data.frame(rbind(Acovctrl, freqAcovctrl))
	
	Acovctrlres = data.frame(rbind(Acovctrl))
	colnames(Acovctrlres) = c("bias2", "var", "mse", "cov1", "cov2")
	print(Acovctrlres)

	res = Acovctrlres
}


resfilename = paste(c(estctrl, trial, randseed, dataseed, factorseed, "res.csv"), 
	collapse="_")
write.table(res, paste(c(savedir, resfilename), 
	collapse="/"), sep = ",", qmethod = "double")

