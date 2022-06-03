
//multiple linear regression code
data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
  //int<lower=0> mu_beta; // prior mu for beta
  //int<lower=0> sig_beta; // prior sig for beta
  real mu_beta[K]; // prior mu for beta
  real sig_beta[K]; // prior sig for beta
  real<lower=0> a; // prior IG parameter - sigma
  real<lower=0> b; // prior IG parameter - sigma
}
parameters {
  vector[K] beta;       // coefficients for predictors (including intercept, if applicable)
  real<lower=0> sigma;  // error scale
}
model {
  // Observational model
  y ~ normal(x * beta, sigma);  // likelihood
  
  // Prior model
  for(j in 1:K){
    beta[j] ~ normal(mu_beta[j], sig_beta[j]);
  }
  sigma ~ inv_gamma(a, b); 
}

