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
from datetime import *
import argparse
import random

import utils
reload(utils)
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('-nc', '--numcauses', \
    type=int, default=100000)
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
    type=int, default=2000)
parser.add_argument('-fm', '--factormodel', \
    choices=['PCA','PPCA','PF','LFA','GMM', 'DEF'], default="DEF")

args, unknown = parser.parse_known_args()

n_causes = args.numcauses
n_units = args.numunits
simset = args.simset
n_iter = args.niter
dataseed = args.dataseed
randseed = args.factorseed
factormodel = args.factormodel
alpha = args.alpha


outdir = wfactorfitdir(simset, n_units, n_causes, args.alpha, factormodel)
if not os.path.exists(outdir):
    os.makedirs(outdir)

#############################################################
# set random seed
#############################################################

# randseed = int(time.time()*1000000%100000000)
print("random seed: ", randseed)
random.seed(randseed)
np.random.seed(randseed)
tf.set_random_seed(randseed)

#############################################################
# load simulated genes
#############################################################

simdatdir = wsimdatdir(simset, n_units, n_causes, args.alpha)
G = np.load(simdatdir+'/'+str(dataseed)+'_snps.npy')
lambdas = np.load(simdatdir+'/'+str(dataseed)+'_groups.npy')

# holdout data
x_train, x_vad, holdout_mask = holdout_data(G)

if factormodel == "PCA":
    # the stochastic vi code subsamples on columns. we pass in the
    # transpose of x_train to subsampling on rows.

    pca = PCA(n_components=50)
    pca.fit(x_train)  
    pca_z_post_np = pca.fit_transform(G)
    pca_x_post_np = pca.inverse_transform(pca_z_post_np)   

    print(pca_z_post_np.shape)
    print(pca_x_post_np.shape)

    print("check pca fit")
    print("trivial mse", ((G-0)**2).mean())
    print("pca mse", ((G-pca_x_post_np)**2).mean())

    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', pca_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', pca_x_post_np)

elif factormodel == "PPCA":

    # the stochastic vi code subsamples on columns. we pass in the
    # transpose of x_train to subsampling on rows.

    ppca_x_post, ppca_w_post, ppca_z_post, \
        ppca_x_post_np, ppca_z_post_np = \
        fit_ppca(x_train.T, 1-holdout_mask.T, \
            stddv_datapoints=1.0, M=100, K=50, \
            n_iter=n_iter, optimizer="adam")

    print("check PPCA fit")

    print("trivial mse", ((G-0)**2).mean())

    print("PPCA mse", ((G-ppca_x_post_np.T)**2).mean())

    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', ppca_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', ppca_x_post_np.T)

    ppca_pval = ppca_predictive_check(x_train.T, x_vad.T, holdout_mask.T, ppca_x_post, ppca_w_post, ppca_z_post)

    print("PPCA predictive check", ppca_pval)
    np.savetxt(outdir+'/'+str(randseed)+'_pval.csv', np.array(ppca_pval))

elif factormodel == "PF":
    # the stochastic vi code subsamples on columns. we pass in the
    # transpose of x_train to subsampling on rows.

    pmf_x_post, pmf_z_post, pmf_w_post, \
        pmf_x_post_np, pmf_z_post_np = \
        fit_pmf(x_train.T, 1-holdout_mask.T, \
            M=100, K=50, n_iter=n_iter, optimizer="adam")    

    print("check PMF fit")

    print("trivial mse", ((G-0)**2).mean())

    print("PMF mse", ((G-pmf_x_post_np.T)**2).mean())

    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', pmf_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', pmf_x_post_np.T)

    pmf_pval = pmf_predictive_check(x_train.T, x_vad.T, holdout_mask.T, pmf_x_post, pmf_w_post, pmf_z_post)

    print("PMF predictive check", pmf_pval)
    np.savetxt(outdir+'/'+str(randseed)+'_pval.csv', np.array(pmf_pval))

elif factormodel == "DEF":

    # the stochastic vi code subsamples on columns. we pass in the
    # transpose of x_train to subsampling on rows. also need to return
    # a different def_z_post_np in utils.py. please see explanation in
    # utils.py in the fit_def() function.

    # in training def, please make sure the loss is negative,
    # otherwise it is often thanks to optimization failure that makes
    # the learning of def fail and cannot deconfound.

    def_x_post, def_z3_post, def_z2_post, def_z1_post, def_W0_post, \
        def_x_post_np, def_z_post_np = \
        fit_def(x_train.T, 1-holdout_mask.T, K=[100,30,5], \
            prior_a=0.1, prior_b=0.3, \
            optimizer=tf.train.RMSPropOptimizer(1e-3),
            # optimizer="adam", \
            n_iter=n_iter)

    # optimizer=tf.train.AdamOptimizer(1e-2), 

    print("check DEF fit")

    print("trivial mse", ((G-0)**2).mean())

    print("DEF mse", ((G-def_x_post_np.T)**2).mean())

    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', def_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', def_x_post_np.T)

    def_pval = def_predictive_check(x_train.T, x_vad.T, holdout_mask.T, def_x_post, def_z1_post, def_W0_post)

    print("DEF predictive check", def_pval)  
    np.savetxt(outdir+'/'+str(randseed)+'_pval.csv', np.array(def_pval))


elif factormodel == "LFA":
    # the stochastic vi code subsamples on columns. we pass in the
    # transpose of x_train to subsampling on rows.

    lfa_x_post, lfa_z_post, lfa_w_post, \
        lfa_x_post_np, lfa_z_post_np = \
        fit_lfa(x_train.T, 1-holdout_mask.T, \
            M=100, K=50, n_iter=n_iter)

    print("check LFA fit")

    print("trivial mse", ((G-0)**2).mean())

    print("LFA mse", ((G-lfa_x_post_np.T)**2).mean())

    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', lfa_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', lfa_x_post_np.T)

    lfa_pval = lfa_predictive_check(x_train.T, x_vad.T, holdout_mask.T, lfa_x_post, lfa_w_post, lfa_z_post)

    print("LFA predictive check", lfa_pval)
    np.savetxt(outdir+'/'+str(randseed)+'_pval.csv', np.array(lfa_pval))


elif factormodel == "GMM":
    gmm_x_post, gmm_mu_post, gmm_sigmasq_post, gmm_pi_post, \
        gmm_x_post_np, gmm_z_post_np = \
        fit_gmm(x_train.T, 1-holdout_mask.T, \
            M=100, K=10, n_iter=n_iter)

    print("check GMM fit")

    print("trivial mse", ((G-0)**2).mean())

    print("GMM mse", ((G-gmm_x_post_np.T)**2).mean())


    np.save(outdir+'/'+str(randseed)+'_Zhat.npy', gmm_z_post_np)
    np.save(outdir+'/'+str(randseed)+'_Ahat.npy', gmm_x_post_np.T)

    gmm_pval = gmm_predictive_check(x_train.T, x_vad.T, holdout_mask.T, gmm_x_post, gmm_mu_post, gmm_sigmasq_post, gmm_pi_post, K=10)

    print("GMM predictive check", gmm_pval)
    np.savetxt(outdir+'/'+str(randseed)+'_pval.csv', np.array(gmm_pval))