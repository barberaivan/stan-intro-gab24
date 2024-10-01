data {
  int<lower=0> n;
  array[n] int y;
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  // Previa
  target += uniform_lpdf(theta | 0, 1);

  // Likelihood
  target += bernoulli_lpmf(y | theta);
}