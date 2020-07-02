# python 2
# edward 1.3.5

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tensorflow as tf
import edward as ed
import numpy as np
import numpy.random as npr
import pandas as pd
import math
import os
import itertools
import time
from decimal import Decimal
import errno
from sklearn import linear_model
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from edward.models import Normal, Gamma, Dirichlet, InverseGamma, \
    Poisson, PointMass, Empirical, ParamMixture, \
    MultivariateNormalDiag, Categorical, Laplace, \
    MultivariateNormalTriL, Bernoulli, TransformedDistribution,\
    Binomial, MultivariateNormalFullCovariance
from edward.util import Progbar
from scipy import sparse, stats
from scipy.special import expit, logit, logsumexp
from sklearn.metrics import r2_score, accuracy_score
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def next_batch(x_train, M):
    # subsample M columns
    D, N = x_train.shape
    idx_batch = np.random.choice(N, M)
    return x_train[:, idx_batch], idx_batch

def wsimdatdir(simset, n_units, n_causes, alpha):
    if simset in ["BN", "TGP", "HGDP"]:
        simdatdir = os.path.join("../dat/simdat", simset, 'units'+str(n_units), 'causes'+str(n_causes))
    elif simset in ["PSD", "SP"]:
        simdatdir = os.path.join("../dat/simdat", simset, 'alpha'+str(alpha), 'units'+str(n_units), 'causes'+str(n_causes))
    return simdatdir

def wfactorfitdir(simset, n_units, n_causes, alpha, factormodel):
    if not os.path.exists("../res"):
        try:
            os.makedirs("../res", 0o700)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if not os.path.exists("../res/factorfit"):
        try:
            os.makedirs("../res/factorfit", 0o700)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    if simset in ["BN", "TGP", "HGDP"]:
        factorfitdir = os.path.join("../res/factorfit", simset, 'units'+str(n_units), 'causes'+str(n_causes), 'model'+factormodel)
    elif simset in ["PSD", "SP"]:
        factorfitdir = os.path.join("../res/factorfit", simset, 'alpha'+str(alpha), 'units'+str(n_units), 'causes'+str(n_causes), 'model'+factormodel)
    return factorfitdir


def woutfile(simset, n_units, n_causes, alpha, randseed, n_iter, cv, snpsig, confint, causalprop, aslin, aslog, aalin, aalog):
    if not os.path.exists("../res/causalest"):
        try:
            os.makedirs("../res/causalest", 0o700)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise    
    if simset in ["BN", "TGP", "HGDP"]:
        causalestfile = "../res/causalest/"+ simset + '_units'+str(n_units) + '_causes'+str(n_causes) + "_randseed" + str(randseed) + "_nitr" + str(n_iter) + "_cv" + str(cv) + "_snpsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "_aslin" + str(format_e(Decimal(aslin))) + "_aslog" + str(format_e(Decimal(aslog))) + "_aalin" + str(format_e(Decimal(aalin))) + "_aalog" + str(format_e(Decimal(aalog)))
    elif simset in ["PSD", "SP"]:
        causalestfile = "../res/causalest/"+ simset + 'alpha'+str(alpha) + '_units'+str(n_units) + '_causes'+str(n_causes) + "_randseed" + str(randseed) + "_nitr" + str(n_iter) + "_cv" + str(cv) + "_snpsig" + str(snpsig) + "_confint" + str(confint) + "_cp" + str(causalprop) + "_aslin" + str(format_e(Decimal(aslin))) + "_aslog" + str(format_e(Decimal(aslog))) + "_aalin" + str(format_e(Decimal(aalin))) + "_aalog" + str(format_e(Decimal(aalog)))
    return causalestfile


def sim_genes_BN(Fs, ps, n_causes, n_units, n_genes, D=3):
    idx = npr.randint(n_genes, size = n_causes)
    p = ps[idx]
    F = Fs[idx]
    Gammamat = np.zeros((n_causes, D))
    for i in range(D):
        Gammamat[:,i] = npr.beta((1-F)*p/F, (1-F)*(1-p)/F)
    S = npr.multinomial(1, (60/210, 60/210, 90/210), size = n_units)
    F = S.dot(Gammamat.T)
    G = npr.binomial(2, F)
    scaledS = G[:,npr.randint(G.shape[1], size=10)]
    lambdas = KMeans(n_clusters=3, random_state=42).fit(scaledS).labels_
    return G, lambdas

def sim_genes_TGP(n_causes, n_units, tgp_pc, D=3):
    S = np.column_stack([tgp_pc])
    scaler = MinMaxScaler()
    scaler.fit(S)
    S = scaler.transform(S)
    Gammamat = np.zeros((n_causes, 3))
    Gammamat[:,0] = 0.45*npr.uniform(size=n_causes)
    Gammamat[:,1] = 0.45*npr.uniform(size=n_causes)
    Gammamat[:,2] = 0.05*np.ones(n_causes)
    S = np.column_stack((S[npr.choice(S.shape[0],size=n_units,replace=True),], \
        np.ones(n_units)))
    F = S.dot(Gammamat.T)
    G = npr.binomial(2, F)
    scaledS = G[:,npr.randint(G.shape[1], size=10)]
    lambdas = KMeans(n_clusters=3, random_state=123).fit(scaledS).labels_
    return G, lambdas

def sim_genes_HGDP(n_causes, n_units, tgp_pc, D=3):
    return sim_genes_TGP(n_causes, n_units, tgp_pc, D=3)

def sim_genes_PSD(Fs, ps, n_causes, n_units, n_genes, D=3, alpha=0.5):
    idx = npr.randint(n_genes, size = n_causes)
    p = ps[idx]
    F = Fs[idx]
    Gammamat = np.zeros((n_causes, D))
    for i in range(D):
        Gammamat[:,i] = npr.beta((1-F)*p/F, (1-F)*(1-p)/F)
    S = npr.dirichlet((alpha, alpha, alpha), size = n_units)
    F = S.dot(Gammamat.T)
    G = npr.binomial(2, F)
    scaledS = G[:,npr.randint(G.shape[1], size=10)]
    lambdas = KMeans(n_clusters=3, random_state=123).fit(scaledS).labels_
    # sG = sparse.csr_matrix(G)
    return G, lambdas

def sim_genes_SP(Fs, ps, n_causes, n_units, D=3, a=0.1):
    Gammamat = np.zeros((n_causes, 3))
    Gammamat[:,0] = 0.45*npr.uniform(size=n_causes)
    Gammamat[:,1] = 0.45*npr.uniform(size=n_causes)
    Gammamat[:,2] = 0.05*np.ones(n_causes)
    S = npr.beta(a, a, size=(n_units, 2))
    S = np.column_stack((S, np.ones(n_units)))
    F = S.dot(Gammamat.T)
    G = npr.binomial(2, F)
    scaledS = G[:,npr.randint(G.shape[1], size=10)]
    lambdas = KMeans(n_clusters=3, random_state=123).fit(scaledS).labels_
    return G, lambdas
    

def sim_single_traits(lambdas, G, a, b, causalprop=0.05, base_scale=10, bin_scale=0.5):
    tau = 1./npr.gamma(3, 1, size=3)
    sigmasqs = tau[lambdas]
    epsilons = npr.normal(0, sigmasqs)
    n_causes = G.shape[1]
    betas = npr.normal(0, 0.5**2, size=n_causes)
    causal_snps = int(causalprop*n_causes)
    betas[causal_snps:] = 0.0

    c = 1 - a - b

    raw_y = G.dot(betas)
    raw_y_std = raw_y.std()

    y = raw_y + \
        np.sqrt(b)*raw_y_std/np.sqrt(a)*(lambdas-lambdas.mean())/lambdas.std() + \
        np.sqrt(c)*raw_y_std/np.sqrt(a)*epsilons/epsilons.std()

    scale = raw_y_std / base_scale
    y = (y - y.mean()) / scale

    y_bin = npr.binomial(1, expit(y / bin_scale))

    true_lambdas = \
        np.sqrt(b)*raw_y_std/np.sqrt(a)*(lambdas)/lambdas.std() / scale
    true_betas = betas / scale

    print("confounding strength np.corrcoef(y, true_lambdas)", \
        np.corrcoef(y, true_lambdas)[0,1])

    return y, y_bin, true_betas, true_lambdas

def fit_outcome_linear(X, y, true_betas, n_causes, alpha=0, CV=False, verbose=False):
    if CV == True:
        linear_reg = regression_CV(X, y, outtype="linear", verbose=verbose)
    else:
        linear_reg = regression_noCV(X, y, alpha, outtype="linear", verbose=verbose)
    
    linear_rmse = np.sqrt(((true_betas - linear_reg.coef_[:n_causes])**2).mean())
    trivial_rmse = np.sqrt(((true_betas - 0)**2).mean())
    if verbose:
        print("linear outcome rmse", linear_rmse, "\nlinear - trivial", linear_rmse - trivial_rmse)
    return linear_reg, linear_rmse


def fit_outcome_logistic(X, y_bin, true_betas, n_causes, alpha=0, CV=False, verbose=False):
    if CV == True:
        logistic_reg = regression_CV(X, y_bin, outtype="logistic", verbose=verbose)
    else:
        logistic_reg = regression_noCV(X, y_bin, alpha, outtype="logistic", verbose=verbose)     
    logistic_rmse = np.sqrt(((true_betas - logistic_reg.coef_[0][:n_causes])**2).mean())
    trivial_rmse = np.sqrt(((true_betas - 0)**2).mean())
    if verbose:
        print("logistic outcome rmse: ", logistic_rmse, \
            "\nlogistic - trivial: ", logistic_rmse - trivial_rmse)
    return logistic_reg, logistic_rmse


def holdout_data(X, holdout_portion = 0.2):
    # randomly holdout some entries of X
    num_datapoints, data_dim = X.shape

    n_holdout = int(holdout_portion * num_datapoints * data_dim)

    holdout_row = np.random.randint(num_datapoints, size=n_holdout)
    holdout_col = np.random.randint(data_dim, size=n_holdout)
    holdout_mask = (sparse.coo_matrix((np.ones(n_holdout), \
        (holdout_row, holdout_col)), shape = X.shape)).toarray()
    # exclude repeated sampling of holdout entries.
    holdout_mask = np.minimum(holdout_mask, np.ones(X.shape)) 
    holdout_mask = np.float32(holdout_mask)

    holdout_subjects = np.unique(holdout_row)

    x_train = np.multiply(1-holdout_mask, X)
    x_vad = np.multiply(holdout_mask, X)
    return x_train, x_vad, holdout_mask


def fit_ppca(x_train, holdout_mask, stddv_datapoints=1.0, M=100, K=10, n_iter=5000, optimizer=tf.train.RMSPropOptimizer(1e-4)):


    # the following code subsets on the column of x_train
    D, N = x_train.shape

    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    # MODEL
    idx_ph = tf.placeholder(tf.int32, M)
    x_ph = tf.placeholder(tf.float32, [D, M])
    holdout_ph = tf.placeholder(tf.float32, [D, M])

    w = Normal(loc=0.0, scale=10.0, sample_shape=[D, K])
    z = Normal(loc=0.0, scale=1.0, sample_shape=[M, K])
    x = Normal(loc=tf.multiply(tf.matmul(w, z, transpose_b=True), holdout_ph),
               scale=stddv_datapoints*tf.ones([D, M]))

    # INFERENCE
    qw_variables = [tf.Variable(tf.random_normal([D, K])),
                    tf.Variable(tf.random_normal([D, K]))]

    qw = Normal(loc=qw_variables[0], scale=tf.nn.softplus(qw_variables[1]))

    qz_variables = [tf.Variable(tf.random_normal([N, K])),
                    tf.Variable(tf.random_normal([N, K]))]

    qz = Normal(loc=tf.gather(qz_variables[0], idx_ph),
                scale=tf.nn.softplus(tf.gather(qz_variables[1], idx_ph)))

    inference_w = ed.KLqp({w: qw}, data={x: x_ph, z: qz})
    inference_z = ed.KLqp({z: qz}, data={x: x_ph, w: qw})

    scale_factor = float(N) / M
    inference_w.initialize(scale={x: scale_factor, z: scale_factor},
                           var_list=qz_variables,
                           n_samples=5, 
                           n_iter=n_iter, optimizer=optimizer)
    inference_z.initialize(scale={x: scale_factor, z: scale_factor},
                           var_list=qw_variables,
                           n_samples=5, optimizer=optimizer)

    sess = ed.get_session()
    tf.global_variables_initializer().run()
    loss = []
    for _ in range(inference_w.n_iter):
        x_batch, idx_batch = next_batch(x_train, M)
        holdout_batch = holdout_mask[:,idx_batch]
        for _ in range(5):
            inference_z.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})

        info_dict = inference_w.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
        inference_w.print_progress(info_dict)

        t = info_dict['t']
        loss.append(info_dict['loss'])

    w_post = Normal(loc=qw_variables[0], scale=tf.nn.softplus(qw_variables[1]))
    z_post = Normal(loc=qz_variables[0],
                scale=tf.nn.softplus(qz_variables[1]))
    x_post = Normal(loc=tf.matmul(w_post, z_post, transpose_b=True),
               scale=stddv_datapoints*tf.ones([D, N]))

    ppca_x_post_np = x_post.mean().eval()
    ppca_z_post_np = z_post.mean().eval()

    return x_post, w_post, z_post, ppca_x_post_np, ppca_z_post_np


def ppca_predictive_check(x_train, x_vad, holdout_mask, x_post, w_post, z_post, stddv_datapoints=1.0, n_rep=10, n_eval=10):
    '''
    n_rep: the number of replicated datasets we generate
    n_eval: the number of samples we draw samples from the inferred Z and W
    '''
    holdout_row, holdout_col = np.where(holdout_mask > 0)
    holdout_gen = np.zeros([n_rep,x_train.shape[0], x_train.shape[1]])

    for i in range(n_rep):
        x_generated = x_post.sample().eval()

        # look only at the heldout entries
        holdout_gen[i] = np.multiply(x_generated, holdout_mask)

    obs_ll = []
    rep_ll = []
    for j in range(n_eval):
        print("eval", j)
        w_sample = w_post.sample().eval()
        z_sample = z_post.sample().eval()
        
        holdoutmean_sample = np.multiply(w_sample.dot(z_sample.T), holdout_mask)
        obs_ll.append(np.mean(stats.norm(holdoutmean_sample, \
                            stddv_datapoints).logpdf(x_vad), axis=0))

        rep_ll.append(np.mean(stats.norm(holdoutmean_sample, \
                            stddv_datapoints).logpdf(holdout_gen),axis=1))
        
    obs_ll_per_zi, rep_ll_per_zi = np.mean(np.array(obs_ll), axis=0), np.mean(np.array(rep_ll), axis=0)

    pvals = np.array([np.mean(rep_ll_per_zi[:,i] < obs_ll_per_zi[i]) for i in range(len(obs_ll_per_zi))])
    holdout_subjects = np.unique(holdout_col)
    overall_pval = np.mean(pvals[holdout_subjects])
    print("Predictive check p-values", overall_pval)
    return overall_pval


def fit_pmf(x_train, holdout_mask, gamma_prior=0.1, M=100, K=10, n_iter=20000, optimizer=tf.train.RMSPropOptimizer(1e-4)):
    # here the holdout mask is an indicator of the training entries.

    # the following code subsets on the column of x_train
    D, N = x_train.shape

    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    idx_ph = tf.placeholder(tf.int32, M)
    x_ph = tf.placeholder(tf.float32, [D, M])
    holdout_ph = tf.placeholder(tf.float32, [D, M])

    U = Gamma(gamma_prior, gamma_prior, sample_shape=[M,K])
    V = Gamma(gamma_prior, gamma_prior, sample_shape=[D,K])
    x = Poisson(tf.matmul(V, U, transpose_b=True))

    min_scale = 1e-5

    qV_variables = [tf.Variable(tf.random_uniform([D, K])), \
                    tf.Variable(tf.random_uniform([D, K]))]

    qV = TransformedDistribution(
                distribution=Normal(qV_variables[0],\
                                    tf.maximum(tf.nn.softplus(qV_variables[1]), \
                                               min_scale)),
                bijector=tf.contrib.distributions.bijectors.Exp())


    qU_variables = [tf.Variable(tf.random_uniform([N, K])), \
                    tf.Variable(tf.random_uniform([N, K]))]


    qU = TransformedDistribution(\
                distribution=Normal(tf.gather(qU_variables[0], idx_ph),\
                    tf.maximum(tf.nn.softplus(tf.gather(qU_variables[1], idx_ph)), min_scale)),
                bijector=tf.contrib.distributions.bijectors.Exp())


    scale_factor = float(N) / M

    # We apply variational EM with E-step over local variables
    # and M-step to point estimate the global weight matrices.
    inference_e = ed.KLqp({U: qU},
                        data={x: x_ph, V:qV})
    inference_m = ed.KLqp({V:qV},
                       data={x: x_ph, U:qU})

    inference_e.initialize(scale={x: scale_factor, U: scale_factor}, var_list=qU_variables, optimizer=optimizer)
    inference_m.initialize(scale={x: scale_factor, U: scale_factor}, optimizer=optimizer)

    tf.global_variables_initializer().run()

    loss = []
    n_iter_per_epoch = 1000
    n_epoch = int(n_iter / n_iter_per_epoch)
    for epoch in range(n_epoch):
        print("Epoch {}".format(epoch))
        nll = 0.0

        pbar = Progbar(n_iter_per_epoch)
        for t in range(n_iter_per_epoch):
            x_batch, idx_batch = next_batch(x_train, M)
            holdout_batch = holdout_mask[:,idx_batch]
    #         x_batch = x_batch.todense()
            pbar.update(t)
            info_dict_e = inference_e.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            info_dict_m = inference_m.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            nll += info_dict_e['loss'] 
            loss.append(info_dict_e['loss'])
        print('epoch', epoch, '\n loss: ', str(info_dict_e['loss']))   

    V_post = TransformedDistribution(
            distribution=Normal(qV_variables[0],\
                                tf.maximum(tf.nn.softplus(qV_variables[1]), \
                                           min_scale)),
            bijector=tf.contrib.distributions.bijectors.Exp())
    U_post = TransformedDistribution(
                distribution=Normal(qU_variables[0],\
                                    tf.maximum(tf.nn.softplus(qU_variables[1]), \
                                               min_scale)),
                bijector=tf.contrib.distributions.bijectors.Exp())
    x_post = Poisson(tf.matmul(V_post, U_post, transpose_b=True))

    pmf_x_post_np = x_post.mean().eval()
    pmf_z_post_np = U_post.eval()

    return x_post, U_post, V_post, pmf_x_post_np, pmf_z_post_np


def pmf_predictive_check(x_train, x_vad, holdout_mask, x_post, V_post, U_post,n_rep=10, n_eval=10):
    '''
    n_rep: the number of replicated datasets we generate
    n_eval: the number of samples we draw samples from the inferred Z and W
    '''
    holdout_row, holdout_col = np.where(holdout_mask > 0)
    holdout_gen = np.zeros([n_rep, x_train.shape[0], x_train.shape[1]])

    for i in range(n_rep):
        x_generated = x_post.sample().eval()

        # look only at the heldout entries
        holdout_gen[i] = np.multiply(x_generated, holdout_mask)

    obs_ll = []
    rep_ll = []
    for j in range(n_eval):
        U_sample = U_post.sample().eval()
        V_sample = V_post.sample().eval()
        
        holdoutmean_sample = np.multiply(V_sample.dot(U_sample.T), holdout_mask)
        obs_ll.append(\
            np.mean(np.ma.masked_invalid(stats.poisson.logpmf(np.array(x_vad, dtype=int), \
                holdoutmean_sample)), axis=0))

        rep_ll.append(\
            np.mean(np.ma.masked_invalid(stats.poisson.logpmf(holdout_gen, \
                holdoutmean_sample)), axis=1))
        
    obs_ll_per_zi, rep_ll_per_zi = np.mean(np.array(obs_ll), axis=0), np.mean(np.array(rep_ll), axis=0)

    pvals = np.array([np.mean(rep_ll_per_zi[:,i] < obs_ll_per_zi[i]) for i in range(len(obs_ll_per_zi))])
    holdout_subjects = np.unique(holdout_col)
    overall_pval = np.mean(pvals[holdout_subjects])
    print("Predictive check p-values", overall_pval)

    return overall_pval



def fit_def(x_train, holdout_mask, K=[100,30,5], M=100, prior_a=0.1, prior_b=0.3, shape=0.1, q='lognormal', optimizer=tf.train.RMSPropOptimizer(1e-4), n_iter=20000):


    # here the holdout mask is an indicator of the training entries.

    # we default to RMSProp here. but we can also use Adam (lr=1e-2).
    # A successful training of def usually means a negative ELBO. In
    # this code, we used the stopping criterion being two consecutive
    # iterations of positive ELBO change. Alternative, we can stop at
    # the iteration that is 10% larger than the minimum ELBO.

    # this code subsample on row. if we want to subsample on columns,
    # then we can pass in a transpose of the original matrix. all else
    # stay the same. and return transpose of def_x_post_np and also let
    # def_z_post_np = qz3_post.eval().T (check the dimensionality to
    # make sure the number of rows == the number of units.)

    # the same trick applies to all other factor models.

    # subsample on rows
    D, N = x_train.shape # number of documents, vocabulary size
    logdir = '~/log/def/'
    logdir = os.path.expanduser(logdir)
    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    idx_ph = tf.placeholder(tf.int32, M)
    x_ph = tf.placeholder(tf.float32, [D, M])
    holdout_ph = tf.placeholder(tf.float32, [D, M])
    
    # MODEL
    W2 = Gamma(prior_a, prior_b, sample_shape=[K[2], K[1]])
    W1 = Gamma(prior_a, prior_b, sample_shape=[K[1], K[0]])
    W0 = Gamma(prior_a, prior_b, sample_shape=[K[0], M])

    z3 = Gamma(prior_a, prior_b, sample_shape=[D, K[2]])
    z2 = Gamma(shape, shape / tf.matmul(z3, W2))
    z1 = Gamma(shape, shape / tf.matmul(z2, W1))
    x = Poisson(tf.matmul(z1, W0))


    # INFERENCE
    def pointmass_q(shape, tfvar=None):
        min_mean = 1e-3
        mean_init = tf.random_normal(shape)
        if tfvar is None:
            rv = PointMass(tf.maximum(tf.nn.softplus(tf.Variable(mean_init)), min_mean))
        else:
            mean = tfvar
            rv = PointMass(tf.maximum(tf.nn.softplus(mean), min_mean))
        return rv


    def gamma_q(shape):
        # Parameterize Gamma q's via shape and scale, with softplus constraints.
        min_shape = 1e-3
        min_scale = 1e-5
        shape_init = 0.5 + 0.1 * tf.random_normal(shape)
        scale_init = 0.1 * tf.random_normal(shape)
        rv = Gamma(tf.maximum(tf.nn.softplus(tf.Variable(shape_init)),
                            min_shape),
                 tf.maximum(1.0 / tf.nn.softplus(tf.Variable(scale_init)),
                            1.0 / min_scale))
        return rv


    def lognormal_q(shape, tfvar=None):
        min_scale = 1e-5
        loc_init = tf.random_normal(shape)
        scale_init = 0.1 * tf.random_normal(shape)
        if tfvar is None:
            rv = TransformedDistribution(
              distribution=Normal(
                  tf.Variable(loc_init),
                  tf.maximum(tf.nn.softplus(tf.Variable(scale_init)), min_scale)),
              bijector=tf.contrib.distributions.bijectors.Exp())
        else:
            loctfvar, scaletfvar = tfvar
            rv = TransformedDistribution(
              distribution=Normal(
                  loctfvar,
                  tf.maximum(tf.nn.softplus(scaletfvar), min_scale)),
              bijector=tf.contrib.distributions.bijectors.Exp())
        return rv

    qW0all = tf.Variable(tf.random_normal([K[0], N]))
    qW0sub = tf.gather(qW0all, idx_ph, axis=1)

    qW2 = pointmass_q(W2.shape)
    qW1 = pointmass_q(W1.shape)
    qW0 = pointmass_q(W0.shape, qW0sub)
    if q == 'gamma':
        qz3 = gamma_q(z3.shape)
        qz2 = gamma_q(z2.shape)
        qz1 = gamma_q(z1.shape)
    else:
        qz3 = lognormal_q(z3.shape)
        qz2 = lognormal_q(z2.shape)
        qz1 = lognormal_q(z1.shape)

    # We apply variational EM with E-step over local variables
    # and M-step to point estimate the global weight matrices.
    inference_e = ed.KLqp({z1: qz1, z2: qz2, z3: qz3},
                          data={x: x_ph, W0: qW0, W1: qW1, W2: qW2})
    inference_m = ed.MAP({W0: qW0, W1: qW1, W2: qW2},
                         data={x: x_ph, z1: qz1, z2: qz2, z3: qz3})

    timestamp = str(int(time.time()*1e7))
    logdir += timestamp + '_' + '_'.join([str(ks) for ks in K]) + \
        '_q_' + str(q)
    kwargs = {'optimizer': optimizer,
              'n_print': 100,
              'logdir': logdir,
              'log_timestamp': False}
    if q == 'gamma':
        kwargs['n_samples'] = 30
    inference_e.initialize(**kwargs)
    inference_m.initialize(optimizer=optimizer)

    tf.global_variables_initializer().run()


    n_iter_per_epoch = 1000
    n_epoch = int(n_iter / n_iter_per_epoch)
    print("n_epoch", n_epoch)
    min_nll = 1e16
    prev_change = -1e16
    for epoch in range(n_epoch):
        print("Epoch {}".format(epoch))
        nll = 0.0

        pbar = Progbar(n_iter_per_epoch)
        for t in range(1, n_iter_per_epoch + 1):
            x_batch, idx_batch = next_batch(x_train, M)
            holdout_batch = holdout_mask[:,idx_batch]
            pbar.update(t)
            info_dict_e = inference_e.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            info_dict_m = inference_m.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            nll += info_dict_e['loss']

        # Compute perplexity averaged over a number of training iterations.
        # The model's negative log-likelihood of data is upper bounded by
        # the variational objective.
        nll = nll / n_iter_per_epoch
        perplexity = np.exp(nll / np.sum(x_train))
        print("Negative log-likelihood <= {:0.3f}".format(nll))
        print("Perplexity <= {:0.3f}".format(perplexity))
        z3_post = qz3
        z2_post = Gamma(shape, shape / tf.matmul(z3_post, W2))
        z1_post = Gamma(shape, shape / tf.matmul(z2_post, W1))
        W0_post = pointmass_q(qW0all.shape, qW0all)
        x_post = Poisson(tf.matmul(z1_post, W0_post))
        def_x_post_np = x_post.mean().eval()
        print("trivial mse", ((x_train-0)**2).mean())
        print("mse", ((x_train-def_x_post_np)**2).mean())
        print(nll, min_nll, nll < min_nll)
        if nll < min_nll:
            min_nll = nll.copy()
            min_z3_post = z3_post
            min_z2_post = z2_post
            min_z1_post = z1_post
            min_W0_post = W0_post
            min_x_post = x_post
            min_def_x_post_np = def_x_post_np.copy()
            min_def_z_post_np = W0_post.eval().T.copy()

        cur_change = (nll - min_nll)/np.abs(min_nll)

        print("cur-LL", nll, "min-LL", min_nll, "diffratio-LL", cur_change)

        print(prev_change, cur_change)

        if prev_change > 0:
            if cur_change > 0:
                break

        prev_change = cur_change

        if math.isnan(nll):
            break

    z3_post = qz3
    z2_post = Gamma(shape, shape / tf.matmul(z3_post, W2))
    z1_post = Gamma(shape, shape / tf.matmul(z2_post, W1))
    W0_post = pointmass_q(qW0all.shape, qW0all)
    x_post = Poisson(tf.matmul(z1_post, W0_post))
    def_x_post_np = x_post.mean().eval()
    def_z_post_np = W0_post.eval().T

    return x_post, z3_post, z2_post, z1_post, W0_post, min_def_x_post_np, min_def_z_post_np



def def_predictive_check(x_train, x_vad, holdout_mask, x_post, z1_post, W0_post, n_rep=10, n_eval=10):
    '''
    n_rep: the number of replicated datasets we generate
    n_eval: the number of samples we draw samples from the inferred Z and W
    '''
    holdout_row, holdout_col = np.where(holdout_mask > 0)
    holdout_gen = np.zeros([n_rep, x_train.shape[0], x_train.shape[1]])

    for i in range(n_rep):
        x_generated = x_post.sample().eval()

        # look only at the heldout entries
        holdout_gen[i] = np.multiply(x_generated, holdout_mask)

    obs_ll = []
    rep_ll = []
    for j in range(n_eval):
        z1_sample = z1_post.sample().eval()
        W0_sample = W0_post.eval()
        
        holdoutmean_sample = np.multiply(z1_sample.dot(W0_sample), holdout_mask)
        obs_ll.append(\
            np.mean(np.ma.masked_invalid(stats.poisson.logpmf(np.array(x_vad, dtype=int), \
                holdoutmean_sample)), axis=0))

        rep_ll.append(\
            np.mean(np.ma.masked_invalid(stats.poisson.logpmf(holdout_gen, \
                holdoutmean_sample)), axis=1))
        
    obs_ll_per_zi, rep_ll_per_zi = np.mean(np.array(obs_ll), axis=0), np.mean(np.array(rep_ll), axis=0)

    pvals = np.array([np.mean(rep_ll_per_zi[:,i] < obs_ll_per_zi[i]) for i in range(len(obs_ll_per_zi))])
    holdout_subjects = np.unique(holdout_col)
    overall_pval = np.mean(pvals[holdout_subjects])
    print("Predictive check p-values", overall_pval)

    return overall_pval


def fit_lfa(x_train, holdout_mask, M=100, K=10, n_iter=20000, optimizer=tf.train.RMSPropOptimizer(1e-4)):

    # here the holdout mask is an indicator of the training entries.

    D, N = x_train.shape

    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    idx_ph = tf.placeholder(tf.int32, M)
    x_ph = tf.placeholder(tf.float32, [D, M])
    holdout_ph = tf.placeholder(tf.float32, [D, M])

    U = Normal(loc=tf.zeros([K, D]), scale=tf.ones([K, D]))
    V = Normal(loc=tf.zeros([K, M]), scale=tf.ones([K, M]))
    x = Binomial(total_count=2 * tf.ones([D, M]), \
             probs=tf.sigmoid(tf.matmul(U, V, transpose_a=True)), \
             value=tf.zeros([D, M], dtype=tf.float32))
    
    # INFERENCE

    qV_variable = tf.Variable(0.01*tf.random_normal([N, K]))


    # INFERENCE
    qU = PointMass(tf.Variable(tf.random_normal([K, D])))
    qV = PointMass(tf.transpose(tf.gather(qV_variable, idx_ph)))

    inference_e = ed.MAP({U: qU}, data={x: x_ph, V: qV})
    inference_m = ed.MAP({V: qV}, data={x: x_ph, U: qU})
    
    scale_factor = float(N)/M
    optimizer_e = optimizer
    optimizer_m = optimizer
    inference_e.initialize(scale={x: scale_factor, U: scale_factor}, optimizer=optimizer_e)
    inference_m.initialize(scale={x: scale_factor, U: scale_factor}, optimizer=optimizer_m)


    tf.global_variables_initializer().run()

    loss = []
    n_iter_per_epoch = 1000
    n_epoch = int(n_iter / n_iter_per_epoch)

    for epoch in range(n_epoch):
        print("Epoch {}".format(epoch))
        pbar = Progbar(n_iter_per_epoch)
        for t in range(n_iter_per_epoch):
            x_batch, idx_batch = next_batch(x_train, M)
            holdout_batch = holdout_mask[:,idx_batch]
            pbar.update(t)
            info_dict_e = inference_e.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            info_dict_m = inference_m.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch, holdout_ph: holdout_batch})
            loss.append(info_dict_e['loss'])

        print('epoch', epoch, '\n loss: ', str(info_dict_e['loss']))   

    U_post = qU
    V_post = Normal(loc=tf.transpose(qV_variable), scale=tf.ones([K,N]))
    x_post = Binomial(total_count=2 * tf.ones([D, N]), \
                 probs=tf.sigmoid(tf.matmul(U_post, V_post, transpose_a=True)), \
                 value=tf.zeros([D, N], dtype=tf.float32))

    lfa_x_post_np = x_post.mean().eval()
    lfa_z_post_np = qV_variable.eval()

    return x_post, U_post, V_post, lfa_x_post_np, lfa_z_post_np




def lfa_predictive_check(x_train, x_vad, holdout_mask, x_post, V_post, U_post,n_rep=10, n_eval=10):
    '''
    n_rep: the number of replicated datasets we generate
    n_eval: the number of samples we draw samples from the inferred Z and W
    # the column of x_train matrix is subject
    '''
    holdout_row, holdout_col = np.where(holdout_mask > 0)
    holdout_gen = np.zeros([n_rep, x_train.shape[0], x_train.shape[1]])


    # ---- sampling method for Binomial ---
    def _sample_n(self, n=1, seed=None):
        '''
        modified from http://edwardlib.org/api/model-development
        '''
        np_sample = lambda n, p, size: np.random.binomial(n=n, p=p, size=size).astype(np.int32)    
        
        # assume total_count and probs have same as shape
        val = tf.py_func(np_sample, 
                         [tf.cast(self.total_count, tf.int32), self.probs, self.probs.shape], 
                         [tf.int32])[0]
        
        batch_event_shape = self.batch_shape.concatenate(self.event_shape)
        shape = tf.concat([tf.expand_dims(n, 0), tf.convert_to_tensor(batch_event_shape)], 0)
        val = tf.reshape(val, shape)
        
        return val

    Binomial._sample_n = _sample_n

    for i in range(n_rep):
        x_generated = x_post.sample().eval()

        # look only at the heldout entries
        holdout_gen[i] = np.multiply(x_generated, holdout_mask)

    obs_ll = []
    rep_ll = []
    for j in range(n_eval):
        U_sample = U_post.sample().eval()
        V_sample = V_post.sample().eval()
        
        holdoutmean_sample = np.multiply(expit(U_post.sample().eval().T.dot(V_post.sample().eval())), holdout_mask)

        obs_ll.append(np.mean(np.ma.masked_invalid(stats.binom.logpmf(np.array(x_vad, dtype=int), 2, holdoutmean_sample)), axis=0))

        rep_ll.append(np.mean(np.ma.masked_invalid(stats.binom.logpmf(holdout_gen, 2, holdoutmean_sample)), axis=1))
        
    obs_ll_per_zi, rep_ll_per_zi = np.mean(np.array(obs_ll), axis=0), np.mean(np.array(rep_ll), axis=0)

    pvals = np.array([np.mean(rep_ll_per_zi[:,i] < obs_ll_per_zi[i]) for i in range(len(obs_ll_per_zi))])
    holdout_subjects = np.unique(holdout_col)
    overall_pval = np.mean(pvals[holdout_subjects])
    print("Predictive check p-values", overall_pval)

    return overall_pval


def fit_gmm(x_train, holdout_mask, M=100, K=10, n_iter=1000):
    # use Gibbs sampling; no subsampling; marginalize out the categorical latent
    # hand compute the categorical latent in the end
    D, N = x_train.shape

    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    # MODEL
    x_ph = tf.placeholder(tf.float32, [M, D])

    pi = Dirichlet(tf.ones(K))
    mu = Normal(tf.zeros(D), tf.ones(D), sample_shape=K)
    sigmasq = InverseGamma(tf.ones(D), tf.ones(D), sample_shape=K)
    x = ParamMixture(pi, {'loc': mu, 'scale_diag': tf.sqrt(sigmasq)}, MultivariateNormalDiag, sample_shape=M)
    z = x.cat

    T = n_iter  # number of MCMC samples
    # T = 500 # to avoid gmm training forever (GMM struggles on high-dimensional data)
    qpi = Empirical(tf.get_variable(
        "qpi/params", [T, K],
        initializer=tf.constant_initializer(1.0 / K)))
    qmu = Empirical(tf.get_variable(
        "qmu/params", [T, K, D],
        initializer=tf.zeros_initializer()))
    qsigmasq = Empirical(tf.get_variable(
        "qsigmasq/params", [T, K, D],
        initializer=tf.ones_initializer()))
    qz = Empirical(tf.get_variable(
        "qz/params", [T, M],
        initializer=tf.zeros_initializer(),
        dtype=tf.int32))

    inference = ed.Gibbs({pi: qpi, mu: qmu, sigmasq: qsigmasq, z: qz},
                     data={x: x_ph})
    inference.initialize()


    tf.global_variables_initializer().run()
    t_ph = tf.placeholder(tf.int32, [])
    running_cluster_means = tf.reduce_mean(qmu.params[:t_ph], 0)

    try:
        inference.update()
    except Exception:
        pass

    for _ in range(inference.n_iter):
        x_batch, idx_batch = next_batch(x_train, M)
        info_dict = inference.update(feed_dict={x_ph: x_batch.T})
        inference.print_progress(info_dict)
        t = info_dict['t']
        if t % inference.n_print == 0:
            print("\nInferred cluster means:")
            print(sess.run(running_cluster_means, {t_ph: t - 1}))

    mu_s = qmu.sample(50).eval()
    sigmasq_s = qsigmasq.sample(500).eval()
    pi_s = qpi.sample(50).eval()

    z_post = np.array([np.log(pi_s.mean(0)[i]) - np.sum((x_train.T-mu_s.mean(0)[i])**2 / sigmasq_s.mean(0)[i], axis=1)
       for i in range(K)]).argmax(axis=0)
    x_post = mu_s.mean(0)[z_post].T
    mu_post = qmu
    sigmasq_post = qsigmasq
    pi_post = qpi

    gmm_x_post_np = mu_s.mean(0)[z_post].T
    gmm_z_post_np = z_post

    return x_post, mu_post, sigmasq_post, pi_post, gmm_x_post_np, gmm_z_post_np


def gmm_predictive_check(x_train, x_vad, holdout_mask, x_post, mu_post, sigmasq_post, pi_post, K=10, n_rep=10, n_eval=10):
    '''
    n_rep: the number of replicated datasets we generate
    n_eval: the number of samples we draw samples from the inferred Z and W
    '''
    holdout_row, holdout_col = np.where(holdout_mask > 0)
    holdout_gen = np.zeros([n_rep,x_train.T.shape[0], x_train.T.shape[1]])

    for i in range(n_rep):
        z_generated = npr.multinomial(n=1, pvals=pi_post.sample(50).eval().mean(0), size=x_train.T.shape[0]).argmax(axis=1)
        x_generated = npr.normal(loc=mu_post.sample(50).eval().mean(0)[z_generated], scale=np.sqrt(sigmasq_post.sample(50).eval().mean(0))[z_generated])

        # look only at the heldout entries
        holdout_gen[i] = np.multiply(x_generated, holdout_mask.T)

    obs_ll = []
    rep_ll = []
    for j in range(n_eval):
        mu_sample = mu_post.sample().eval()
        sigmasq_sample = sigmasq_post.sample().eval()
        pi_sample = pi_post.sample().eval()
        
        obs_ll_cur = logsumexp(np.array([np.log(pi_sample[i]) - np.sum(np.multiply((x_vad.T-mu_sample[i]), holdout_mask.T)**2 / sigmasq_sample[i], axis=1) for i in range(K)]), axis=0)

        rep_ll_cur = np.array([logsumexp(np.array([np.log(pi_sample[i]) - np.sum(np.multiply((holdout_gen[j]-mu_sample[i]), holdout_mask.T)**2 / sigmasq_sample[i], axis=1) for i in range(K)]), axis=0) for j in range(holdout_gen.shape[0])])

        obs_ll.append(obs_ll_cur)

        rep_ll.append(rep_ll_cur)
        
    obs_ll_per_zi, rep_ll_per_zi = np.mean(np.array(obs_ll), axis=0), np.mean(np.array(rep_ll), axis=0)

    pvals = np.array([np.mean(rep_ll_per_zi[:,i] < obs_ll_per_zi[i]) for i in range(len(obs_ll_per_zi))])
    holdout_subjects = np.unique(holdout_col)
    overall_pval = np.mean(pvals[holdout_subjects])
    print("Predictive check p-values", overall_pval)

    return overall_pval


def regression_CV(X, y, outtype="linear", verbose=False):

    X_train, X_test, y_train, y_test = \
            train_test_split(X, y, test_size=0.2, random_state=123)

    if outtype == "linear":
        reg = linear_model.Ridge(normalize=True)
    elif outtype == "logistic":
        reg = linear_model.RidgeClassifier(normalize=True)

    parameters = dict(alpha=np.logspace(-5, 5, 10))
    clf = GridSearchCV(reg, parameters)
    clf.fit(X_train, y_train)

    bestalpha = clf.best_estimator_.alpha
    print("best alpha", bestalpha)

    if outtype == "linear":
        reg = linear_model.Ridge(alpha=bestalpha, normalize=True)
    elif outtype == "logistic":
        reg = linear_model.RidgeClassifier(alpha=bestalpha, normalize=True)

    reg.fit(X_train, y_train)
    
    if verbose:
        print("training score", reg.score(X_train, y_train))
        print("predictive score", reg.score(X_test, y_test))

    return reg


def regression_noCV(X, y, alpha=0, outtype="linear", verbose=False):

    X_train, X_test, y_train, y_test = \
            train_test_split(X, y, test_size=0.2, random_state=123)

    if outtype == "linear":
        reg = linear_model.Ridge(alpha=alpha, normalize=True)
    elif outtype == "logistic":
        reg = linear_model.RidgeClassifier(alpha=alpha, normalize=True)

    reg.fit(X_train, y_train)
    
    if verbose:
        print("training score", reg.score(X_train, y_train))
        print("predictive score", reg.score(X_test, y_test))

    return reg

def fit_lmm(X, y, cov, outtype="linear", M=100, n_iter=1000, optimizer="adam", verbose=False):
    # following EMMAX to compute the kinship matrix
    # Psi_ij = 1/S \sum_{s=1}^S \frac{(G_i^s-2 * p_s)(G_j^s-2 * p_s)}{2 p_s (1-p_s)}
    # p_s is the sample average allele frequency at location s
    # G_i^s is the allele frequency for individual i and location s
    # Psi is the sample covariance matrix of the SNP matrix

    # ps = np.array([X[:,s].mean()/2. for s in range(X.shape[1])])

    # covs = np.array([(X[:,s] - 2*p)[:,np.newaxis].dot((X[:,s] - 2*p)[:,np.newaxis].T)/(2*p*(1-p)) for s, p in enumerate(ps) if p!=0]).mean(axis=0)

    # cov = np.cov(X)

    # X_train, X_vad, Y_train, Y_vad = \
    #         train_test_split(X, y, test_size=0.2, random_state=123)
    X_train = X.copy()
    Y_train = y.copy()

    N = X_train.shape[0]
    D = X_train.shape[1]
    assert N == Y_train.shape[0]

    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    x_ph = tf.placeholder(tf.float32, [M, D])
    idx_ph = tf.placeholder(tf.int32, M)
    cov_ph = tf.placeholder(tf.float32, [M, M])

    if outtype=="logistic":
        y_ph = tf.placeholder(tf.int32, [M])
    else:
        y_ph = tf.placeholder(tf.float32, [M])

    b = Normal(loc=tf.zeros([1]), scale=10.*tf.ones([1]))
    w = Normal(loc=tf.zeros([D, 1]), scale=0.1*tf.ones([D, 1]))

    ymean = tf.matmul(x_ph, w) + b 

    sigma_sp = Normal(loc=tf.zeros([1]), scale=10.*tf.ones([1]))
    sigma_dp = Normal(loc=tf.zeros([1]), scale=10.*tf.ones([1]))
    sigma_s = tf.sqrt(tf.exp(sigma_sp))
    sigma_d = tf.sqrt(tf.exp(sigma_dp))

    y = MultivariateNormalFullCovariance(loc=ymean[:,0], covariance_matrix=sigma_s*tf.eye(M) + sigma_d*cov_ph)

    if outtype=="logistic":
        y = Bernoulli(logits=y)
    else:
        y = y

    qsigma_s_vars = [tf.Variable(-1.*tf.random_normal([1]))]
    qsigma_s = PointMass(params=qsigma_s_vars[0])
    qsigma_d_vars = [tf.Variable(-1.*tf.random_normal([1]))]
    qsigma_d = PointMass(params=qsigma_d_vars[0])

    qw_variables = [tf.Variable(-1.*tf.random_normal([D, 1]))]
    qw = PointMass(params=qw_variables[0])

    qb_variables = [tf.Variable(-1.*tf.random_normal([1]))]
    qb = PointMass(params=qb_variables[0])

    all_variables = list(itertools.chain(qw_variables, qb_variables, qsigma_s_vars, qsigma_d_vars))

    scale_factor = 1. * N / M 

    inference = ed.MAP({w: qw, b: qb, sigma_sp: qsigma_s, sigma_dp: qsigma_d}, data={y: y_ph})

    inference.initialize(\
        scale={y: scale_factor}, \
        var_list=all_variables, optimizer=optimizer)

    tf.global_variables_initializer().run()

    loss = []

    for i in range(n_iter):
        idx_batch = np.random.choice(N, M)
        x_batch = X_train[idx_batch]
        y_batch = Y_train[idx_batch]
        cov_batch = cov[:,idx_batch][idx_batch,:]
        info_dict = inference.update(\
            feed_dict={x_ph: x_batch, y_ph: y_batch, cov_ph: cov_batch, idx_ph: idx_batch})
        loss.append(info_dict["loss"])
        if verbose:
            if i % 1000 == 0:
                print("Itr", i, "loss", info_dict["loss"])
        if math.isnan(info_dict["loss"]):
            break

        if verbose:
            if (i + 0) % 1000 == 0:

                wval = qw_variables[0].eval()
                bval = qb_variables[0].eval()

                print("#####################")
                print("Itr", i, "training accuracy")

                y_trainpredmean = bval + np.dot(X_train, wval)
                evalacc(Y_train[:,np.newaxis], y_trainpredmean, outtype)


    wval = qw_variables[0].eval()
    bval = qb_variables[0].eval()

    loss = np.array(loss)
    if verbose:
        print("#####################")
        print("final ",  "training accuracy")

        y_trainpredmean = bval + np.dot(X_train, wval)
        trainacc = evalacc(Y_train[:,np.newaxis], y_trainpredmean, outtype)

    return wval

def evalacc(Y_vad, yvadpredmean, outtype):
    if outtype=="logistic":
        yvadpredmean = expit(yvadpredmean)
        yvadpred = npr.binomial(n=1, p=yvadpredmean)
    else:
        yvadpred = yvadpredmean

    yvadpred, yvadpredmean = np.squeeze(np.array(yvadpred)), np.squeeze(np.array(yvadpredmean))


    if outtype !="logistic":
        print("trivial accuracy", r2_score(Y_vad, np.array([np.ones_like(Y_vad[:,j])*np.mean(Y_vad[:,j]) for j in range(Y_vad.shape[1])]).T))
        print("realized accuracy (hand)", r2_score(Y_vad, yvadpredmean))
        vadacc = r2_score(Y_vad, yvadpredmean)
    else:
        print("trivial accuracy", accuracy_score(Y_vad, np.array([np.ones_like(Y_vad[:,j])*stats.mode(Y_vad[:,j])[0][0] for j in range(Y_vad.shape[1])]).T))
        print("realized accuracy (hand)", accuracy_score(Y_vad, yvadpred)) 
        vadacc = accuracy_score(Y_vad, yvadpred)

    return vadacc

