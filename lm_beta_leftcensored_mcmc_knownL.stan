
//multiple linear regression code
data {
  int<lower=0> N_obs;           // Number of observed data points
  int<lower=0> N_cens;          // Number of censored data points
  int<lower=0> K;               // Number of predictors
  matrix[N_obs, K] X_obs;       // Observed predictor matrix
  matrix[N_cens, K] X_cens;     // Censored predictor matrix
  real y_obs[N_obs];            // Observations (non-censored)
  real<upper=min(y_obs)> L;     // Value where censoring occurs
  real mu_beta[K];              // prior mu for beta
  real sig_beta[K];             // prior sig for beta
  real<lower=0> a;              // prior IG parameter - sigma
  real<lower=0> b;              // prior IG parameter - sigma
}
parameters {
  vector[K] beta;               // coefficients for predictors (including intercept, if applicable)
  real<lower=0> sigma;          // error scale
  real<upper=L> y_cens[N_cens]; // The underlying censored values
}
model {
  // model
  y_obs ~ normal(X_obs * beta, sigma);   // target density - observed points
  y_cens ~ normal(X_cens * beta, sigma); // target density - censored points
  
  // Prior model
  for(j in 1:K){
    beta[j] ~ normal(mu_beta[j], sig_beta[j]);
  }
  sigma ~ inv_gamma(a, b); 
}

