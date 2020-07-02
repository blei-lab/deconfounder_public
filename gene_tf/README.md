# GWAS Simulation

## data

We include a demo dataset in `dat/rawdat/hapmap/`, which is a subset
of the HapMap dataset from
https://github.com/mimno/admixture-ppc/tree/master/hapmap3-files

## environment

python 2.7

tensorflow 1.5.0

edward 1.3.5

## scripts

`src/run_script_simdat.sh` simulates the SNPs and the outcomes.

`src/run_script_fitfactor.sh` fits factor models.

`src/run_script_causalest.sh` estimates causal coefficients.

`src/agg_dfs.py` aggregates the simulation results.

## output

`linear_rmse` and `logistic_rmse` is the RMSE for regression with all
the SNPs together. (lower is better.)

`linear_rmse_sing` and `logistic_rmse_sing` is the RMSE for regression
with one SNP at a time. (lower is better.)

## on fitting factor models with variational Bayes

To fit a factor model, we need to fit it multiple times with different
initialization and select the fit with the highest evidence lower
bound (ELBO). The reason is that variational Bayes performs a
non-convex optimization to fit factor models; it can get stuck in
different local optima with different initializations. Operationally,
it means to fit factor models multiple times with different random
seeds and select the fit with the highest ELBO.


