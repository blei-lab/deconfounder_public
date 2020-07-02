# python 2
# edward 1.3.5

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tensorflow as tf
import edward as ed
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import numpy.random as npr
import pandas as pd
import math
import os
from datetime import *
from sklearn import linear_model
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from edward.models import Normal, Gamma, Dirichlet, InverseGamma, \
    Poisson, PointMass, Empirical, ParamMixture, \
    MultivariateNormalDiag, Categorical, Laplace, \
    MultivariateNormalTriL, Bernoulli, TransformedDistribution,\
    Binomial
from edward.util import Progbar
from scipy import sparse, stats
from scipy.special import expit, logit
from sklearn.metrics import r2_score, accuracy_score
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split




def next_batch(x_train, M):
    # subsample M columns
    D, N = x_train.shape
    idx_batch = np.random.choice(N, M)
    return x_train[:, idx_batch], idx_batch


def holdout_data(X):
    # randomly holdout some entries of X
    num_datapoints, data_dim = X.shape

    holdout_portion = 0.2
    n_holdout = int(holdout_portion * num_datapoints * data_dim)

    holdout_row = np.random.randint(num_datapoints, size=n_holdout)
    holdout_col = np.random.randint(data_dim, size=n_holdout)
    holdout_mask = (sparse.coo_matrix((np.ones(n_holdout),                             (holdout_row, holdout_col)),                             shape = X.shape)).toarray()
    holdout_mask = np.minimum(holdout_mask, np.ones(X.shape))
    holdout_mask = np.float32(holdout_mask)


    holdout_subjects = np.unique(holdout_row)

    x_train = np.multiply(1-holdout_mask, X)
    x_vad = np.multiply(holdout_mask, X)
    return x_train, x_vad, holdout_mask


def fit_def(x_train, K=[100,30,5], M=100, prior_a=0.1, prior_b=0.3, shape=0.1, q='lognormal', optimizer=tf.train.RMSPropOptimizer(1e-4), n_iter=20000):

    # we default to RMSProp here. but we can also use Adam (lr=1e-2).
    # A successful training of def usually means a negative ELBO. In
    # this code, we used the stopping criterion being two consecutive
    # iterations of positive ELBO change. Alternative, we can stop at
    # the iteration that is 10% larger than the minimum ELBO.

    # this code subsample on row. if we want to subsample on columns,
    # then we can pass in a transpose of the original matrix. all else
    # stay the same. and return tranpose of def_x_post_np and also let
    # def_z_post_np = qz3_post.eval().T (check the dimensionality to
    # make sure the number of rows == the number of units.)

    # the same trick applies to all other factor models.

    # subsample on rows
    N, D = x_train.shape # number of documents, vocabulary size
    logdir = '~/log/def/'
    logdir = os.path.expanduser(logdir)
    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    idx_ph = tf.placeholder(tf.int32, M)
    x_ph = tf.placeholder(tf.float32, [N, M])
    
    # MODEL
    W2 = Gamma(prior_a, prior_b, sample_shape=[K[2], K[1]])
    W1 = Gamma(prior_a, prior_b, sample_shape=[K[1], K[0]])
    W0 = Gamma(prior_a, prior_b, sample_shape=[K[0], M])

    z3 = Gamma(prior_a, prior_b, sample_shape=[N, K[2]])
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
        # Parameterize Gamma q's via shape and scale, with softplus unconstraints.
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


    # qz3loc, qz3scale = tf.Variable(tf.random_normal([N, K[2]])), tf.Variable(tf.random_normal([N, K[2]]))
    # qz3locsub, qz3scalesub = tf.gather(qz3loc, idx_ph), tf.gather(qz3scale, idx_ph) 

    qW0all = tf.Variable(tf.random_normal([K[0], D]))
    qW0sub = tf.gather(qW0all, idx_ph, axis=1)

    qW2 = pointmass_q(W2.shape)
    qW1 = pointmass_q(W1.shape)
    qW0 = pointmass_q(W0.shape, qW0sub)
    if q == 'gamma':
        # qz3 = gamma_q(z3.shape, (qz3locsub, qz3scalesub))
        qz3 = gamma_q(z3.shape)
        qz2 = gamma_q(z2.shape)
        qz1 = gamma_q(z1.shape)
    else:
        # qz3 = lognormal_q(z3.shape, (qz3locsub, qz3scalesub))
        qz3 = lognormal_q(z3.shape)
        qz2 = lognormal_q(z2.shape)
        qz1 = lognormal_q(z1.shape)

    # We apply variational EM with E-step over local variables
    # and M-step to point estimate the global weight matrices.
    inference_e = ed.KLqp({z1: qz1, z2: qz2, z3: qz3},
                          data={x: x_ph, W0: qW0, W1: qW1, W2: qW2})
    inference_m = ed.MAP({W0: qW0, W1: qW1, W2: qW2},
                         data={x: x_ph, z1: qz1, z2: qz2, z3: qz3})

    timestamp = datetime.strftime(datetime.utcnow(), "%Y%m%d_%H%M%S")
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
    min_nll = 1e16
    prev_change = -1e16
    for epoch in range(n_epoch):
        print("Epoch {}".format(epoch))
        nll = 0.0

        pbar = Progbar(n_iter_per_epoch)
        for t in range(1, n_iter_per_epoch + 1):
            x_batch, idx_batch = next_batch(x_train, M)
            pbar.update(t)
            info_dict_e = inference_e.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch})
            info_dict_m = inference_m.update(feed_dict={x_ph: x_batch, idx_ph: idx_batch})
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
        def_z_post_np = W0_post.eval().T
        print("trivial mse", np.square(x_train).mean())
        print("mse", np.square(x_train-def_x_post_np).mean())
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

        # if cur_change > 0.1:
            # if nll < 0:
            # break
        if math.isnan(nll):
            break


    # z3_post = lognormal_q([N, K[2]], (qz3loc, qz3scale))
    # z2_post = Gamma(shape, shape / tf.matmul(z3_post, W2))
    # z1_post = Gamma(shape, shape / tf.matmul(z2_post, W1))
    # x_post = Poisson(tf.matmul(z1_post, W0))

    # def_x_post_np = x_post.mean().eval()
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
