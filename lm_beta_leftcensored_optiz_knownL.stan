
//multiple linear regression code
data {
  int<lower=0> N_obs;           // Number of observed data points
  int<lower=0> N_cens;          // Number of censored data points
  int<lower=0> K;               // Number of predictors
  matrix[N_obs, K] X_obs;       // Observed predictor matrix
  matrix[N_cens, K] X_cens;     // Censored predictor matrix
  real y_obs[N_obs];            // Observations (non-censored)
  real<upper=min(y_obs)> L;     // Value where censoring occurs
  //real<lower=max(y_obs)> U;    // Value where right censoring occurs
  real mu_beta[K];              // prior mu for beta
  real sig_beta[K];             // prior sig for beta
  real<lower=0> a;              // prior IG parameter - sigma
  real<lower=0> b;              // prior IG parameter - sigma
}
parameters {
  vector[K] beta;               // coefficients for predictors (including intercept, if applicable)
  real<lower=0> sigma;          // error scale
}
model {
  // model
  y_obs ~ normal(X_obs * beta, sigma);   // target density - observed points
  target += N_cens * normal_lcdf(L | X_cens * beta, sigma); // left censored
  //target += N_cens * normal_lccdf(U | X_cens * beta, sigma); // right censored
  
  // Prior model
  for(j in 1:K){
    beta[j] ~ normal(mu_beta[j], sig_beta[j]);
  }
  sigma ~ inv_gamma(a, b); 
}

