rm(list = ls())

rawdir = "./smoking_R"

setwd(paste(c(rawdir, "src"), collapse="/"))

dir.create(paste(c(rawdir, "dat/simdat"), collapse="/"))


###########################################################
# parse arguments
###########################################################


args <- commandArgs(trailingOnly = TRUE)

ntrials = as.integer(args[1])
confscale = as.numeric(args[2])
outcome_model = args[3]
simset = args[4]
randseed = as.integer(args[5])

# ntrials = 30
# confscale = 1
# outcome_model = "linear"
# simset = "indep"
# randseed = 1587247895


if(simset=="indep"){
	savedir = paste(c(rawdir, "dat/simdat/indep"), collapse="/")
	dir.create(savedir)
} else if(simset=="dep"){
	# for causes with dependence	
	savedir = paste(c(rawdir, "dat/simdat/dep"), collapse="/")
	dir.create(savedir)
}

params <- data.frame(
	"ntrials" = c(ntrials), 
	"confscale" = c(confscale), 
	"outcome_model" = c(outcome_model),
	"simset" = c(simset),
	"randseed" = c(randseed))

params_filenames = paste(c(randseed, "simdat_params.csv"), collapse="_")
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
set.seed(randseed)
print(randseed)

###########################################################
# load data
###########################################################

data = read.table("../dat/nmes_data.csv", sep=",", header=TRUE)
num = dim(data)[1]
# this line does not subsample because the total number is <10000;  if
# you want to subsample, change 10000 to be the desired subsample size
subsamp_size = min(10000, num) 
subsamp = sample(1:num, subsamp_size, replace=FALSE)
dat = data[subsamp,]

###########################################################
# set causes and confounder
###########################################################

if(simset=="indep"){
	cau_names = c("packyears","marital")
	cau = dat[cau_names]
} else if(simset=="dep"){
	# for causes with dependence
	cau_names = c("packyears","marital")
	cau = dat[cau_names]
	cauplus = dat["marital"] + rnorm(dim(data)[1])*0.1
	cau = cbind(cau, cauplus)		
}

conf = dat["LASTAGE"]

cov_names = c("AGESMOKE", "MALE", "RACE3", "beltuse", 
	"educate")
cov = dat[cov_names]
cov = data.matrix(cov)


# transform and standardize variables
cau = apply(cau, 2, scale)
conf = scale(conf)
cov = apply(cov, 2, scale) 

# set causes
A = data.matrix(cau)
C = data.matrix(conf)

N = dim(A)[1]
D = dim(A)[2]
Kc = dim(C)[2]

###########################################################
# simulate data
###########################################################

default_val = NA

intercepts = matrix(default_val, ntrials, 1)
betass = matrix(default_val, ntrials, D)
gammass = matrix(default_val, ntrials, Kc)
simmeans = matrix(default_val, ntrials, N)
simYs = matrix(default_val, ntrials, N)

if(outcome_model=="linear"){
	family = gaussian()
} else if (outcome_model=="logistic"){
	family = binomial()
}

causes = A
confounders = C
covariates = cov

for (i in 1:ntrials){
	simdata = sim_conf_data(A, C, outdir, family, 
		confscale=confscale, sd=1.0, 
		seed=randseed) 
	intercepts[i,] = simdata$intercept
	betass[i,] = simdata$betas
	gammass[i,] = simdata$gammas
	simYs[i,] = simdata$simY
	simmeans[i,] = simdata$simmean
}

save_simdat(causes, savedir, randseed, simset)
save_simdat(confounders, savedir, randseed, simset)
save_simdat(covariates, savedir, randseed, simset)
save_simdat(intercepts, savedir, randseed, simset)
save_simdat(betass, savedir, randseed, simset)
save_simdat(gammass, savedir, randseed, simset)
save_simdat(simYs, savedir, randseed, simset)
save_simdat(simmeans, savedir, randseed, simset)

