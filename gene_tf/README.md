# GWAS Simulation

## Data

We include a demo dataset in `dat/rawdat/hapmap/`. It is a subset of
the HapMap dataset from
https://github.com/mimno/admixture-ppc/tree/master/hapmap3-files.

## Environment

python 2.7

tensorflow 1.5.0

edward 1.3.5

## How to execute the scripts

1. Run the script `src/run_script_simdat.sh` to simulate the SNPs and
   the outcomes.

2. Run the script `src/run_script_fitfactor.sh` to fit factor models.

3. With a selected factor model (see the note below), run the script
   `src/run_script_causalest.sh` to estimate causal coefficients.

4. Run the script `src/agg_dfs.py` to aggregate the simulation
   results.


## Output

The files `res/*_avg.csv` contain output from the demo dataset.

The columns `linear_rmse` and `logistic_rmse` in the output files are
the RMSE for regression with all the SNPs together. (Lower is better.)

The columns `linear_rmse_sing` and `logistic_rmse_sing` in the output
files are the RMSE for regression with one SNP at a time. (Lower is
better.)

The file `res/output.md` summarizes the results.

## Differences from Wang and Blei (2019)

The results here differ from Tables 6-15 in Wang and Blei (2019).
However, the captions of those tables apply to these results as well.
One reason for the difference is that we here study a subset of the
data.  There are also several other differences:

1. This implementation fits all factor models with the Edward
   implementation of black box variational inference. The original
   implementation used coordinate accent variational inference (CAVI)
   for Poisson factorization.

2. This implementation uses Adam as the default optimizer. The
   original implementation used RMSprop.

3. This implementation uses the ridge regression in `scikit-learn` to
   estimate causal effects. The original implementation performed a
   Bayesian linear regression with black box variational inference, as
   is implemented in Edward.

## A note on fitting and selecting the factor models

We fit the factor models with variational inference.  To fit a factor
model, we fit it multiple times with different random seeds; each one
leads to a different local optimum of the non-convex optimization of
the variational objective.  In this application, we suggest selecting
the factor model fit that has a high evidence lower bound (ELBO).
