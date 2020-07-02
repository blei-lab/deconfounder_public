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
import errno
from datetime import *
import argparse
import random
import time

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
parser.add_argument('-seed', '--seed', \
    type=int, default=52744889)


args, unknown = parser.parse_known_args()

n_causes = args.numcauses
n_units = args.numunits
simset = args.simset
randseed = args.seed
alpha = args.alpha / 100.


#############################################################
# set random seed
#############################################################

# randseed = int(time.time()*1000000%100000000)
print("random seed: ", randseed)
random.seed(randseed)
np.random.seed(randseed)
tf.set_random_seed(randseed)


if not os.path.exists("../dat"):
    print("No raw data!")

if not os.path.exists("../dat/simdat"):
    try:
        os.makedirs("../dat/simdat", 0o700)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
            
#############################################################
# simulate genes (causes) and traits (outcomes)
#############################################################
# load Hapmap data for BN, PSD, SP
# to preprocess the data, run clean_hapmap.py
Fs = np.loadtxt("../dat/rawdat/hapmap/clean_csv/Fs.csv")
ps = np.loadtxt("../dat/rawdat/hapmap/clean_csv/ps.csv")
genes = pd.read_csv("../dat/rawdat/hapmap/clean_csv/genes.csv")
n_genes = genes.shape[1]

if simset == "BN":
    G, lambdas = sim_genes_BN(Fs, ps, n_causes, n_units, n_genes)
elif simset == "TGP":
    tgp_pc = np.array(pd.read_csv("../dat/rawdat/tgp/tgp_pca.csv"))[:, :2]
    G, lambdas = sim_genes_TGP(n_causes, n_units, tgp_pc)
elif simset == "HGDP":
    hgdp = np.array(pd.read_csv("../dat/rawdat/hgdp/hgdp_subset.csv").T[1:])
    pca = PCA(n_components=2, svd_solver='full')
    pca.fit(hgdp)
    hgdp_pc = pca.transform(hgdp)
    G, lambdas = sim_genes_HGDP(n_causes, n_units, hgdp_pc)
elif simset == "PSD":
    G, lambdas = sim_genes_PSD(Fs, ps, n_causes, n_units, n_genes, alpha=alpha)
elif simset == "SP":
    G, lambdas = sim_genes_SP(Fs, ps, n_causes, n_units, a=alpha)


# save simulated data
simdatdir = wsimdatdir(simset, n_units, n_causes, args.alpha)

if not os.path.exists(simdatdir):
    os.makedirs(simdatdir)

# remove genes that take the same value on all individuals
const_cols = np.where(np.var(G,axis=0)<0.001)[0]
print(const_cols)
if len(const_cols) > 0:
    G = G[:,list(set(range(n_causes))-set(const_cols))]
    n_causes -= len(const_cols)
    
print(G.shape)



np.save(simdatdir+'/'+str(randseed)+'_snps.npy', G)
np.save(simdatdir+'/'+str(randseed)+'_groups.npy', lambdas)

# np.savetxt(simdatdir+'/'+str(randseed)+'_snps.csv', G)
# np.savetxt(simdatdir+'/'+str(randseed)+'_groups.csv', lambdas)

