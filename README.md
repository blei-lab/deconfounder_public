# Reference implementation of the deconfounder

This folder contains the code for:

+ the empirical study about smoking (Section 6.1)

+ the empirical study about genome-wide association studies (GWAS)
  (Section 6.2)

+ the empirical study about movies (Section 6.3)

# Acknowledgments

We thank Justin Grimmer, Dean Knox, and Brandon Stewart for pointing
out concerns with a previous version of the code.

+ Smoking study:

	+ There was a bug with how the code handled heldout data in the
	  Stan file for fitting factor models. This bug has been fixed.

	+ Regarding the definition of bias, variance, and MSE metrics: We
	  now have added an explanation of these metrics in
	  `utils.R`. These metrics are defined as average per-simulation
	  posterior bias/variance/MSE, following Korattikara et
	  al. (2014); Chen et al. (2015); Gustafson (2015).

	+ Regarding the number of posterior samples drawn to compute the
	  Monte Carlo estimate of posterior bias/variance/MSE: We have
	  increased the default number of posterior samples to 10. (This
	  number can be set as a hyperparameter and further increased.)

+ Movie study:

	+ [Grimmer et al.  (2020)](https://www.dropbox.com/s/71m4ncw6s9nkek7/gks.pdf)
	  points out that conditioning on observed covariates can be
	  important for causal estimation. Using the movie data, they
	  illustrate how the deconfounder without conditioning on observed
	  covariates can produce unreasonable causal estimates.

      The exploratory results on actors in Wang and Blei (2019) do not
	  condition on covariates.  The "intervened test set" analyses in
	  the supplement of Wang and Blei (2019) includes such
	  conditioning.

      For users that are interested in exploring different analyses,
	  the causal estimation section of
	  `movie_py/movie-actor-py2-causalest.ipynb` includes a block
	  of code that conditions on observed covariates in addition to
	  the substitute confounder.

# References

Chen, C., Ding, N., & Carin, L. (2015). On the convergence of
stochastic gradient MCMC algorithms with high-order integrators. In
Advances in Neural Information Processing Systems (pp. 2278-2286).


Grimmer, J., Knox, & D., Stewart, B. (2020). Naive regression requires
weaker assumptions than factor models to adjust for multiple cause
confounding
[[link](https://www.dropbox.com/s/71m4ncw6s9nkek7/gks.pdf)]


Gustafson, P. (2015). Bayesian inference for partially identified
models: Exploring the limits of limited data (Vol. 140). CRC Press.


Korattikara, A., Chen, Y., & Welling, M. (2014). Austerity in MCMC
land: Cutting the Metropolis-Hastings budget. In International
Conference on Machine Learning (pp. 181-189).


Wang, Y. and Blei, D.M. (2019) The Blessings of Multiple Causes.
_Journal of American Statistical Association_, 114:528, 1574-1596.
[[link](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2019.1686987?af=R)]






