# implementation in R and stan

Extends model using Poisson-Binomial distribution (of which the binomial distribution is general case) to allow sampling from non-identically distributed sub-populations.

As this results in a non-conjugate prior for the likelihood function, so I use MCMC to sample from the posterior.

Dependencies are  `rstan` and `bayesplot`.
