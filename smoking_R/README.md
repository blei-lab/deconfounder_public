# Smoking simulation

## Data

The NMES dataset is in `dat/nmes_data.csv`, taken from the `causaldrf`
package in R.

## Environment

R version 3.5.3 (2019-03-11)

## How to execute the scripts

1. Run the script `src/run_script_simdat.sh` to simulate the outcomes.

2. Run the script `src/run_script_fitfactor.sh` to fit factor models.

3. With a selected factor model (see the note below), run the script
   `src/run_script_causalest.sh` to estimate causal coefficients.

4. Run the script `src/agg_dfs.py` to aggregate the simulation
   results.

## Output

The files `res/*_avg.csv` include output from this implementation.

The file `res/output.md` is a summary of these output files.

## Differences from Wang and Blei (2019)

This output differs from the results in Table 3 of Wang and Blei
(2019). In this output, controlling for the substitute confounder Zhat
performs similarly as controlling for the reconstructed causes Ahat;
and the capacity of factor models or including additional covariates
does not significantly increase the variance. There are two
differences in this implementation:

1. The original implementation used a hand-coded Stan program for
   Bayesian linear regression. This implementation uses the
   `stan_glm.fit` function with `QR=TRUE`, as implemented in
   `rstanarm` package.

2. The original implementation used a Gamma(0.1, 0.1) prior for the
   standard deviation of the regression error and hierarchical Cauchy
   priors for the regression coefficients and intercepts. This
   implementation uses an exponential(1) prior for the standard
   deviation of the error and Cauchy(0, 10) priors for the regression
   coefficients and intercepts.  (The original priors are not easily
   implemented with the rstanarm package.)


## A note on fitting and selecting the factor models

We fit the factor models with variational inference.  To fit a factor
model, we fit it multiple times with different random seeds; each one
leads to a different local optimum of the non-convex optimization of
the variational objective.

Grimmer et al. 2020 point out that different fits of the factor model
can lead to different causal estimates in this study.  With just two
causes, there can be high variance in the downstream estimates.  To
control this variance, we suggest selecting the factor model fit that
trades off several properties: (a) a high evidence lower bound (ELBO)
(b) non-perfect correlation between Zhat (or Ahat) and A, e.g., around
0.5, and (c) low correlation among the components of Zhat or Ahat (so
that they are nearly independent).  The reason for including aspects
of the correlation structure is to lower the variance of the causal
estimates.
