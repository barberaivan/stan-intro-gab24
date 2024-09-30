data {
  int<lower=0> n;
  array[n] int y;
  vector[n] z;
  
  // parámetros que definen las previas (conocidos, por eso van en data)
  real alpha_mu;
  real alpha_sigma;
  real beta_mu;
  real beta_sigma;
}

parameters {
  real alpha;
  real beta;
}

transformed parameters {
  // theta ya no es un parámetro, sino una cantidad derivada, por eso va
  // en este bloque. Y como theta varía con z, ya no es un escalar (real), sino
  // un vector de largo n.
  vector[n] theta = inv_logit(alpha + beta * z);
  // no hace falta restringir a theta explícitamente con
  // vector<lower=0, upper=1>[n] porque inv_logit asegura eso.
}

model {
  // Previas sobre los parámetros, no sobre theta
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(beta_mu, beta_sigma);

  // Likelihood
  y ~ bernoulli(theta);
}