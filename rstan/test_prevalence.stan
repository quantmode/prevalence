// J.A. Goldfrank 2020, MIT License
//
// The input data
data {
  int m; // number of trials
  int N[m]; // number of positive results in each trial m
  int n[m]; // total sampled population in each trial m
  real alpha[m]; // prior alpha for each trial m (beta-binomial)
  real beta[m]; // prior beta for each trial m (beta-binomial)
  real tp_alpha; //test true positive alpha (beta-binomial)
  real tp_beta; //test true positive beta (beta-binomial)
  real fp_alpha; //test false positive alpha (beta-binomial)
  real fp_beta; //test false positive beta (beta-binomial)
}

// The model parameters
parameters {
  real<lower=0, upper=1> p[m];
  real<lower=0, upper=1> tpp;
  real<lower=0, upper=1> fpp;
}

// Transformed parameters (not used)
transformed parameters {
}

// The model to be estimated.
model {
  p ~ beta(alpha, beta); //probability of having disease - vectorized for each m trial
  tpp ~ beta(tp_alpha, tp_beta); //true positive probability
  fpp ~ beta(fp_alpha, fp_beta); //false positive probability
  for (i in 1:m){
    for (j in 1:N[m]){
      target += log_sum_exp(log(tpp) + log(p[i]),
                log(fpp) + log(1-p[i]));
    }
    for (j in 1:n[m]-N[m]){
      target += log_sum_exp(log(1-tpp) + log(p[i]),
                log(1-fpp) + log(1-p[i]));
    }
  } //binomial-poission, equivalent to binomial for single trial
}
