data {
  int n;
  int y[n];
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  // Previa
  theta ~ uniform(0, 1);

  // Likelihood
  y ~ bernoulli(theta);
}