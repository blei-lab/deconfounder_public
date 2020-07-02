# Smoking simulation

## data

The NMES dataset is at `dat/nmes_data.csv`, taken from the
`causaldrf` package in R.

## environment

R version 3.5.3 (2019-03-11)

## bash scripts to run the code

`src/run_script_simdat.sh` simulates the outcomes.

`src/run_script_fitfactor.sh` fits factor models.

`src/run_script_causalest.sh` estimates causal coefficients.

`src/agg_dfs.py` aggregates the simulation results.

## on fitting factor models with variational Bayes

To fit a factor model, we need to fit it multiple times with different
random seeds; the fits will be different because different seed leads
to different initialization in non-convex optimization and gets to a
different local optima.

We usually select the factor model fit with a high evidence lower
bound (ELBO) and a close to 0.5 correlation between the substitute
confounder Zhat and the causes (if we adjust for the substitute
confounder Zhat). If we adjust for the reconstructed causes Ahat, we
would want a close to 0.5 correlation between the reconstructed causes
Ahat and the causes A. High correlations often blow up the variance of
the estimates and lead to large squared error. The error is often
smaller when the correlation is close to 0.5.


<!-- ## on evaluation metrics

consider the estimate \hat{\beta} as a randomly drawn approximate
posterior sample of the beta parameter. we compute the following bias,
variance, and mean squared error:

bias = E[\hat{\beta}] - \beta

var = Var(\hat{\beta})

mse = E[(\hat{\beta} - \beta)^2]

the expectation is over the approximate posterior.

these definitions follow the bias / variance / mse of posterior
samples (Korattikara et al., 2014; Chen et al., 2015).

Korattikara, A., Chen, Y., & Welling, M. (2014). Austerity in MCMC
land: Cutting the Metropolis-Hastings budget. In International
Conference on Machine Learning (pp. 181-189).

Chen, C., Ding, N., & Carin, L. (2015). On the convergence of
stochastic gradient MCMC algorithms with high-order integrators. In
Advances in Neural Information Processing Systems (pp. 2278-2286). -->
