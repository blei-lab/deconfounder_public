# Movie study

## Data

The directory `dat/tmdb_raw/tmdb-5000-movie-dataset/` contains the raw
TMDB 5000 Dataset, downloaded from
https://www.kaggle.com/tmdb/tmdb-movie-metadata.

## Environment

python 2

tensorflow 1.5.0

edward 1.3.5

## Scripts

Run `src/tmdb_movie_preprocess.ipynb` to preprocess data.

Run `src/movie-actor-py2-subconf.ipynb` to fit factor models.

Run `src/movie-actor-py2-causalest.ipynb` to estimate causal effect.

## Output

The file `res/factorfit/` contains prefit substitute confounders.

The file `src/movie-actor-py2-causalest.ipynb` uses these prefit
substitute confounders to estimate causal effects.

## Differences from Wang and Blei (2019)

The results in `src/movie-actor-py2-causalest.ipynb` differ from the
results in Section 6.3 and Tables 16-18 in Wang and Blei (2019).
However, the captions of those tables apply to these results as well.
There are two differences in this implementation:

1. This implementation uses coordinate ascent variational inference
  (CAVI) for Poisson factorization. The original implementation fitted
  all factor models with the Edward implementation of black box
  variational inference.

2. This implementation uses OLS in `scikit-learn` to estimate causal
  effects. The original implementation performed a Bayesian linear
  regression with black box variational inference, as implemented in
  Edward.

## A note on conditioning on observed covariates

[Grimmer et al.
(2020)](https://www.dropbox.com/s/71m4ncw6s9nkek7/gks.pdf) points out
that conditioning on observed covariates can be important for causal
estimation. Using the movie data, they illustrate how the deconfounder
without conditioning on observed covariates can produce unreasonable
causal estimates.

The exploratory results on actors in Wang and Blei (2019) do not
condition on covariates.  The "intervened test set" analyses in the
supplement of Wang and Blei (2019) includes such conditioning.

For users that are interested in exploring different analyses, the
causal estimation section of
`movie_py/movie-actor-py2-causalest.ipynb` includes a block of code
that conditions on observed covariates in addition to the substitute
confounder.
