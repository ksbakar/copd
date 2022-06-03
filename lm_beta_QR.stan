
//multiple linear regression code
data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
  real mu_beta[K]; // prior mu for beta
  real sig_beta[K]; // prior sig for beta
  real<lower=0> a; // prior IG parameter - sigma
  real<lower=0> b; // prior IG parameter - sigma
}
// this step does some transformations to the data
transformed data {
  // QR decomposition
  matrix[N, K] Q; 
  matrix[K, K] R; 
  matrix[K, K] R_inv;
  // 
  Q = qr_Q(x)[, 1:K] * N;
  R = qr_R(x)[1:K, ] / N;
  R_inv = inverse(R);
}
parameters {
  vector[K] theta;      // coefficients on Q
  real<lower=0> sigma;  // error scale
}
transformed parameters {
  vector[K] beta = R_inv * theta;
}
model {
  // likelihood
  //y ~ normal(Q * theta, sigma);  
  target += normal_lpdf(y | Q * theta, sigma);
  
  // Prior model
  for(j in 1:K){
    //beta[j] ~ normal(mu_beta[j], sig_beta[j]);
    target += normal_lpdf(beta[j] | mu_beta[j], sig_beta[j]);
  }
  sigma ~ inv_gamma(a, b); 
}

