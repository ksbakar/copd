
data {
  int <lower = 0> N; // Defining the number of data points
  int <lower = 0, upper = 1> y [N]; // A variable that describes [1] or [0]
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  real mu_beta[K]; // prior for beta - mu
  real sig_beta[K]; // prior for beta - sig
}
parameters {
  vector[K] beta;   // coefficients for predictors (including intercept, if applicable)
}
model {
  // Observational model
  y ~ bernoulli_logit(x * beta); // likelihood

  // Prior model
  for(j in 1:K){
    beta[j] ~ normal(mu_beta[j], sig_beta[j]);
  }
}
