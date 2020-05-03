## R Version 3.6.3
## J.A. Goldfrank 2020, MIT License
##

library(rstan)
library(bayesplot)
rstan_options(auto_write = TRUE)


n <- 3300 # number of people
N <- 50 # number of positives
m <- length(N)  
alpha <- 0.5 # alpha of beta-binomial prior probability on prevalence
beta <- 0.5 # beta of beta-binomial prior probability on prevalence
tp_alpha <- 103 # alpha of beta-binomial prior probability on true positive rate
tp_beta <- 19 # alpha of beta-binomial prior probability on true positive rate
fp_alpha <- 2 # beta of beta-binomial prior probability on true positive rate
fp_beta <- 399 # beta of beta-binomial prior probability on true positive rate

stan_data <- list(N=N, n=n, m=m, alpha=alpha, beta=beta, tp_alpha=tp_alpha, tp_beta=tp_beta, fp_alpha=fp_alpha, fp_beta=fp_beta)

if (m == 1){
  fit1 <- stan(
    file = "test_prevalence_single.stan",
    data = stan_data,
    chains = 8,
    warmup = 15000,
    iter = 30000,
    thin = 1,
    cores = 8, # should not be higher than your number of cores, but won't break anything
    refresh = 0, # set = 1 to umute
    control=list(adapt_delta = 0.95, max_treedepth = 12)
    )
} else{
  fit1 <- stan(
    file = "test_prevalence.stan",
    data = stan_data,
    chains = 8,
    warmup = 15000,
    iter = 30000,
    thin = 1,
    cores = 8,
    refresh = 0, #set = 1 to umute
    control=list(adapt_delta = 0.95, max_treedepth = 12)
  )
}

print(fit1, digits=4)
png(filename="prevalence_density.png")
mcmc_hist(fit1, pars = c("p"), binwidth=0.0001,freq=TRUE)
dev.off()

png(filename="specificty_density.png")
mcmc_hist(fit1, pars = c("fpp"), binwidth=0.0001,freq=TRUE)
dev.off()
