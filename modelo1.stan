data {
  int<lower=0> n;
  array[n] int y;
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