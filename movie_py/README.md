# A case study about movies

## data

`dat/tmdb_raw/tmdb-5000-movie-dataset/` contains the raw TMDB 5000 Dataset downloaded 
from https://www.kaggle.com/tmdb/tmdb-movie-metadata

## environment

python 2

tensorflow 1.5.0

edward 1.3.5

## scripts

`src/tmdb_movie_preprocess.ipynb` preprocesses the data.

`src/movie-actor-py2-subconf.ipynb` constructs substitute confounders.

`src/movie-actor-py2-causalest.ipynb` estimates causal effects.

## factor model fits
`res/factorfit/` contains prefit substitute confounders;
`src/movie-actor-py2-causalest.ipynb` use these prefit ones to estimate
causal effects.

