# python 2.7.12 (default, Dec  4 2017, 14:50:18) 
# [GCC 5.4.0 20160609]
# tensorflow 1.5.0
# edward 1.3.5

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tensorflow as tf
import edward as ed
import numpy as np
import numpy.random as npr
import pandas as pd
import os
import argparse
import random
from sklearn.metrics import log_loss
import time
from decimal import Decimal

import utils
reload(utils)
from utils import *

#############################################################
# set the scale of simulation
#############################################################

parser = argparse.ArgumentParser()
parser.add_argument('-nc', '--numcauses', \
    type=int, default=5000)
parser.add_argument('-nu', '--numunits', \
    type=int, default=5000)
parser.add_argument('-sim', '--simset', \
    choices=['BN','TGP','HGDP','PSD','SP'], default="BN")
parser.add_argument('-alpha', '--alpha', \
    type=int, default=10)
parser.add_argument('-dataseed', '--dataseed', \
    type=int, default=29200017)
parser.add_argument('-factorseed', '--factorseed', \
    type=int, default=29200422)
parser.add_argument('-nitr', '--niter', \
    type=int, default=5000)
parser.add_argument('-cv', '--cv', \
    type=int, default=0)
parser.add_argument('-snpsig', '--snpsig', \
    type=int, default=40)
parser.add_argument('-confint', '--confint', \
    type=int, default=40)
parser.add_argument('-cp', '--causalprop', \
    type=int, default=10)
parser.add_argument('-aslin', '--alpha_sing_lin', \
    type=float, default=4.)
parser.add_argument('-aslog', '--alpha_sing_log', \
    type=float, default=0.)
parser.add_argument('-aalin', '--alpha_all_lin', \
    type=float, default=0.)
parser.add_argument('-aalog', '--alpha_all_log', \
    type=float, default=1e-4)


args, unknown = parser.parse_known_args()


# n_set_causes is the number of causes in the set up, i.e. 5000,
# 10000; it is different from the actual number of causes because we
# remove SNP columns that are constant

n_set_causes = args.numcauses 
n_units = args.numunits
simset = args.simset
CV = args.cv
a = args.snpsig/100.
b = args.confint/100.
causalprop = args.causalprop/100.
n_iter = args.niter
factorseed = args.factorseed

alpha_sing_lin = args.alpha_sing_lin
alpha_sing_log = args.alpha_sing_log
alpha_all_lin = args.alpha_all_lin
alpha_all_log = args.alpha_all_log


#############################################################
# set random seed
#############################################################

# randseed = 52744889
randseed = int(time.time()*1e7%1e8)
print("random seed: ", randseed)
random.seed(randseed)
np.random.seed(randseed)
tf.set_random_seed(randseed)

# load data

simdatdir = wsimdatdir(simset, n_units, n_set_causes, args.alpha)

G = np.load(simdatdir+'/'+str(args.dataseed)+'_snps.npy')
lambdas = np.load(simdatdir+'/'+str(args.dataseed)+'_groups.npy')

# G = np.loadtxt(simdatdir+'/'+str(args.dataseed)+'_snps.csv', dtype=np.int32)
# lambdas = np.loadtxt(simdatdir+'/'+str(args.dataseed)+'_groups.csv', dtype=np.int32)

n_causes = G.shape[1]
n_units = G.shape[0]

bin_scale = 0.5 # can change it to 0.5 for bin_scale

y, y_bin, true_betas, true_lambdas = sim_single_traits(lambdas, G, \
    a=a, b=b, causalprop=causalprop, bin_scale=bin_scale)
# above: true_lambdas are the precise confounding contribution in the
# generated outcome. it is a scaled version of lambdas.

# set up output file
outfile = woutfile(simset, n_units, n_set_causes, args.alpha, randseed, n_iter, CV, args.snpsig, args.confint, args.causalprop, alpha_sing_lin, alpha_sing_log, alpha_all_lin, alpha_all_log) + '_res.log'

all_adjustments = np.array(['none', 'oracle', 'lmm', 'PCA','PPCA','PF','LFA','GMM','DEF'])
res = pd.DataFrame({"adjustments": all_adjustments, \
    "linear_rmse_sing": -1. * np.ones(len(all_adjustments)), \
    "logistic_rmse_sing": -1. * np.ones(len(all_adjustments)), \
    "linear_rmse": -1. * np.ones(len(all_adjustments)), \
    "logistic_rmse": -1. * np.ones(len(all_adjustments))})

################################
# load all fitted factor models
################################

adjustments = []
z_post_nps = []
x_post_nps = []

################################
# load all fitted factor models
################################
for adjustment in ['PCA','PPCA','PF','LFA','GMM','DEF']:
    print_str = "\n#########################\n" + "loading " + adjustment
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)

    factorfitdir = wfactorfitdir(simset, n_units, n_set_causes, args.alpha, adjustment)
    z_post_file = factorfitdir+'/'+str(factorseed)+'_Zhat.npy'
    if os.path.exists(z_post_file):
        z_post_np = np.load(factorfitdir+'/'+str(factorseed)+'_Zhat.npy')
        x_post_np = np.load(factorfitdir+'/'+str(factorseed)+'_Ahat.npy')
        adjustments.append(adjustment)
        z_post_nps.append(z_post_np)
        x_post_nps.append(x_post_np)
    else:
        print_str = "no fitted "+adjustment+" model found"
        print(print_str)
        with open(outfile, 'a') as f:
            f.write(print_str)

adjustments = np.array(adjustments)

#############################################
# individual adjustment: subset deconfounder
#############################################

print_str = '\n#############################\n' + 'individual effect estimation' + '\n#############################\n'
print(print_str)
with open(outfile, 'a') as f:
    f.write(print_str)


for adjustment in (['none', 'oracle'] + list(adjustments)):
    linear_rmse_sing = np.zeros(n_causes)
    linear_coef_sing = np.zeros(n_causes)
    logistic_rmse_sing = np.zeros(n_causes)
    logistic_coef_sing = np.zeros(n_causes)

    print_str = '\n#############################\n' + adjustment + '\n#############################\n'
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)
    
    for j in range(n_causes):
        if j % 1000 == 0:
            print("Calculating for cause ", j, adjustment)
        if adjustment == 'none':
            X = np.column_stack([G[:,j][:,np.newaxis]])
        elif adjustment == 'oracle':
            X = np.column_stack([G[:,j][:,np.newaxis], lambdas])
        elif adjustment == 'lmm':
            cov = np.cov(G)
            X = np.column_stack([G[:,j][:,np.newaxis]])
            linear_coef_sing[j] = fit_lmm(X, y, cov, outtype="linear", M=5, n_iter=1000, optimizer="adam", verbose=False)
            linear_rmse_sing[j] = np.sqrt(((true_betas[j] - linear_coef_sing[j])**2).mean())
            logistic_coef_sing[j] = fit_lmm(X, y_bin, cov, outtype="logistic", M=5, n_iter=1000, optimizer="adam", verbose=False)
            logistic_rmse_sing[j] = np.sqrt(((true_betas[j] / bin_scale - logistic_coef_sing[j])**2).mean())
            continue
        elif adjustment in ['PCA', 'PPCA', 'GMM']:
            idx = np.where(adjustments == adjustment)[0][0]
            z_post_np = z_post_nps[idx]
            x_post_np = x_post_nps[idx]
            X = np.column_stack([G[:,j][:,np.newaxis], z_post_np])
        elif adjustment in ['PF', 'LFA', 'DEF']:
            idx = np.where(adjustments == adjustment)[0][0]
            z_post_np = z_post_nps[idx]
            x_post_np = x_post_nps[idx]
            X = np.column_stack([G[:,j][:,np.newaxis] - x_post_np[:,j][:,np.newaxis]])    
        
        tmp, linear_rmse_sing[j] = fit_outcome_linear(X, y, true_betas[j], 1, alpha=alpha_sing_lin, CV=CV)
        linear_coef_sing[j] = tmp.coef_[0]

        tmp, logistic_rmse_sing[j] = fit_outcome_logistic(X, y_bin, true_betas[j] / bin_scale, 1, alpha=alpha_sing_log, CV=CV)
        logistic_coef_sing[j] = tmp.coef_[0][0]

    print_str = "\nlinear_rmse_sing.mean() " + str(linear_rmse_sing.mean())
    print_str += "\nlinear_norm_sing.mean() " + str(((linear_coef_sing)**2).mean())
    print_str += "\nlinear_prederr_sing.mean() " + str((((G.dot(linear_coef_sing)-y)**2).mean()))
    print_str += "\n\nlogistic_rmse_sing.mean() " + str(logistic_rmse_sing.mean())
    print_str += "\nlogistic_norm_sing.mean() " + str(((logistic_coef_sing)**2).mean())
    posprob = expit(G.dot(logistic_coef_sing))
    negprob = 1 - posprob
    print_str += "\nlogistic_prederr_sing.mean() " + str(log_loss(y_bin, np.row_stack([posprob, negprob]).T))
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)

    res.iloc[(res["adjustments"] == adjustment).values, 2] = linear_rmse_sing.mean()
    res.iloc[(res["adjustments"] == adjustment).values, 4] = logistic_rmse_sing.mean()


#############################################
# joint adjustment: full deconfounder
#############################################

print_str = '\n#############################\n' + 'joint effect estimation' + '\n#############################\n'
print(print_str)
with open(outfile, 'a') as f:
    f.write(print_str)

for adjustment in (['none', 'oracle', 'lmm'] + list(adjustments)):
    
    print_str = '\n#############################\n' + adjustment + '\n#############################\n'
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)

    if adjustment == 'none':
        X = np.column_stack([G])
    elif adjustment == 'oracle':
        X = np.column_stack([G, lambdas])
    elif adjustment == 'lmm':
        cov = np.cov(G)
        X = np.column_stack([G])
        linear_coef = fit_lmm(X, y, cov, outtype="linear", M=5, n_iter=50000, optimizer="adam", verbose=False)
        linear_rmse = np.sqrt(((true_betas - linear_coef)**2).mean())
        print_str = "linear_rmse " + str(linear_rmse)
        print(print_str)
        with open(outfile, 'a') as f:
            f.write(print_str)

        logistic_coef = fit_lmm(X, y_bin, cov, outtype="logistic", M=5, n_iter=50000, optimizer="adam", verbose=False)
        logistic_rmse = np.sqrt(((true_betas / bin_scale - logistic_coef)**2).mean())
        print_str = "logistic_rmse " + str(logistic_rmse)
        print(print_str)
        with open(outfile, 'a') as f:
            f.write(print_str)
        res.iloc[(res["adjustments"] == adjustment).values, 1] = linear_rmse
        res.iloc[(res["adjustments"] == adjustment).values, 3] = logistic_rmse
        continue
    elif adjustment in ['PCA']:
        print_str = 'no PCA adjustment applicable'
        print(print_str)
        with open(outfile, 'a') as f:
            f.write(print_str)
        continue
    elif adjustment in ['LFA', 'GMM']:
        idx = np.where(adjustments == adjustment)[0][0]
        z_post_np = z_post_nps[idx]
        x_post_np = x_post_nps[idx]
        X = np.column_stack([G, z_post_np])
    elif adjustment in ['PPCA', 'PF', 'DEF']:
        idx = np.where(adjustments == adjustment)[0][0]
        z_post_np = z_post_nps[idx]
        x_post_np = x_post_nps[idx]
        X = np.column_stack([G - x_post_np])   
       

    linear_reg, linear_rmse = fit_outcome_linear(X, y, true_betas, n_causes, alpha=alpha_all_lin, CV=CV, verbose=False)

    print_str = "linear_rmse " + str(linear_rmse)
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)

    logistic_reg, logistic_rmse = fit_outcome_logistic(X, y_bin, true_betas / bin_scale, n_causes, alpha=alpha_all_log, CV=CV, verbose=False)

    print_str = "logistic_rmse " + str(logistic_rmse)
    print(print_str)
    with open(outfile, 'a') as f:
        f.write(print_str)


    res.iloc[(res["adjustments"] == adjustment).values, 1] = linear_rmse
    res.iloc[(res["adjustments"] == adjustment).values, 3] = logistic_rmse

print(res)

res.to_csv(woutfile(simset, n_units, n_set_causes, args.alpha, randseed, n_iter, CV, args.snpsig, args.confint, args.causalprop, alpha_sing_lin, alpha_sing_log, alpha_all_lin, alpha_all_log) + "_res"+ str(int(time.time()*1e7%1e8)) + ".csv")
