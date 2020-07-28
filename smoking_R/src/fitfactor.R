rm(list = ls())


rawdir = "smoking_R"

setwd(paste(c(rawdir, "src"), collapse="/"))

dir.create(paste(c(rawdir, "res/factorfit"), collapse="/"))


# Please see README.md for a note on fitting and selecting the factor
# models.


###########################################################
# parse arguments
###########################################################


args <- commandArgs(trailingOnly = TRUE)

factor_model = args[1]
factor_model_K = as.integer(args[2])
factor_model_dstd = as.numeric(args[3])
factor_model_Zstd = as.numeric(args[4])
simset = args[5]
dataseed = as.integer(args[6])
randseed = as.integer(args[7])

# factor_model = "quadratic"
# factor_model_K = 1
# factor_model_dstd = 0.1
# factor_model_Zstd = 2.0
# simset = "indep"
# dataseed = 825172439
# randseed = 20144133
# print(randseed)


savedir = paste(c(rawdir, "res", "factorfit", simset,
	paste(c(factor_model, factor_model_K), collapse="_")),
	collapse="/")

dir.create(paste(c(rawdir, "res", "factorfit", simset),
	collapse="/"))
dir.create(savedir)

params <- data.frame(
	"factor_model" = c(factor_model), 
	"factor_model_K" = c(factor_model_K), 
	"factor_model_dstd" = c(factor_model_dstd),
	"factor_model_Zstd" = c(factor_model_Zstd),
	"simset" = c(simset),
	"dataseed" = c(dataseed),
	"randseed" = c(randseed))

params_filenames = paste(c(dataseed, randseed, "fitfactor_params.csv"), collapse="_")
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

# randseed = as.integer(as.numeric(Sys.time()))
# randseed = 123
set.seed(randseed)
print(randseed)


###########################################################
# load causes
###########################################################

if(simset=="indep"){
	simdatdir = paste(c(rawdir, "dat/simdat/indep"), collapse="/")
} else if(simset=="dep"){
	# for causes with dependence	
	simdatdir = paste(c(rawdir, "dat/simdat/dep"), collapse="/")
}

A = load_simdat("causes", simdatdir, dataseed, simset)
C = load_simdat("confounders", simdatdir, dataseed, simset)

N = dim(A)[1]
D = dim(A)[2]
Kc = dim(C)[2]

###########################################################
# hold out data
###########################################################

# the p-value for the check will be more accurate when the
# holdout_portion is larger

A_holdout = holdout_data(A, holdout_portion=0.10)
holdout_mask = A_holdout$holdout_mask
A_train = A_holdout$x_train
A_vad = A_holdout$x_vad
holdout_row = A_holdout$holdout_row

###########################################################
# setup factor model
###########################################################

K = factor_model_K
data_std = factor_model_dstd
Z_std = factor_model_Zstd

if(factor_model=="linear"){
	factor_file = 'linearfactor_knownvar.stan'
	model = stan_model(file = 'linearfactor_knownvar.stan')
} else if(factor_model=="quadratic"){
	factor_file = 'quadraticfactor_knownvar.stan'
	model = stan_model(file = 'quadraticfactor_knownvar.stan')
}

factor_data = list(N=N, D=D, K=K, X_train=A_train, 
	data_std=data_std, Z_std=Z_std, X_vad=A_vad, 
	holdout_mask=holdout_mask, X_all=A)


###########################################################
# fit factor model using hmc
###########################################################

# factorfit = stan(file = factor_file, data = factor_data)

###########################################################
# fit factor model using vb
###########################################################

factormap = optimizing(model, data = factor_data, as_vector=FALSE, 
	iter=1000, seed=randseed)
factorinit = factormap$par
factorfit = vb(model, data = factor_data, init=factorinit, 
	iter=10000, eta=0.25, adapt_engaged=0, tol_rel_obj=0.01, 
	output_samples=1000, seed=randseed)

la = extract(factorfit)

# posterior predictive check
idv_pvals = matrix(-1, N, 1)
for (i in holdout_row){
	rep_lp = data.frame(lp=as.matrix(la$rep_lp[,i]))
	vad_lp = data.frame(lp=as.matrix(la$vad_lp[,i]))
	rep_lp$label = "rep"
	vad_lp$label = "vad"
	idv_pvals[i] = mean(rep_lp["lp"] <= colMeans(vad_lp["lp"]))
	# print(c("idv check", i, idv_pvals[i]))
}

print(c("idv check", mean(idv_pvals[holdout_row])))


all_rep_lp = data.frame(lp=as.matrix(apply(la$rep_lp, 1, sum)))
all_vad_lp = data.frame(lp=as.matrix(apply(la$vad_lp, 1, sum)))
all_rep_lp$label = "rep"
all_vad_lp$label = "vad"

pval = mean(all_rep_lp["lp"] < colMeans(all_vad_lp["lp"]))
print(c("overall check", pval))

check_lps = rbind(all_rep_lp, all_vad_lp)
# ggplot(check_lps, aes(lp, fill = label)) + geom_density(alpha = 0.2)

# two fits of the same model might have different pvals
# make sure the pval is not too low, e.g. <0.1.

###########################################################
# extract substitute confounder and reconstructed causes
###########################################################

# tips: check the correlation between real causes and Ahat / Zhat
# empirically, correlation~0.5 usually works the best
# larger correlation will lead to high variance
# smaller correlation means C will not be deconfounded

# if using a linear outcome model downstream

Zhat = la$Z
Zhat_mean = apply(Zhat,c(2,3),mean)
print(cor(A, Zhat_mean))

Ahat = la$X_pred
Ahat_mean = apply(Ahat,c(2,3),mean)
print(cor(A, Ahat_mean))

ps = la$ps # propensity score
ps_mean = apply(la$ps, 2, mean)


output_samples = dim(Zhat)[1]

print(output_samples)

for (i in 1:output_samples){
	print(i)
	filenameZ = paste(c(dataseed, randseed, "Zhat", i, ".csv"), collapse="_")
	filenameA = paste(c(dataseed, randseed, "Ahat", i, ".csv"), collapse="_")
	write.table(Zhat[i,,], paste(c(savedir, filenameZ), collapse="/"), sep = ",", qmethod = "double")
	write.table(Ahat[i,,], paste(c(savedir, filenameA), collapse="/"), sep = ",", qmethod = "double")
}


