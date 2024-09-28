data {
  int nt;
  int nf;
  int y[nt*nf];
  int x[nt, nf];

  // parámetros que definen las previas (conocidos, por eso van en data)
  real alpha_mu;
  real alpha_sigma;
  real beta_mu;
  real beta_sigma;
}

parameters {
  real alpha;
  real beta;
  real<lower=0, upper=1> iota;
  real<lower=0, upper=1> delta;
  row_vector<lower=0, upper=1>[nf] z1;
}

transformed parameters {
  matrix[nt, nf] theta;
  matrix[nt, nf] z;     // z[2:nf, ] es una cantidad derivada, el único verdadero
                        // parámetro es z[1, ] = z1

  // formas vectorizadas
  vector[nt*nf] theta_vec;
  vector[nt*nf] z_vec;

  // llenamos matrices z y theta
  z[1, ] = z1;
  theta[1, ] = inv_logit(alpha + beta * z1);

  for(f in 1:nf) {
    for(t in 2:nt) {
      // actualizamos z en base al ataque en el intervalo anterior
      z[t, f] =
        z[t-1, f] +
        iota * (1 - z[t-1, f]) * x[t-1, f] -   // ataque (aumenta)
        delta * z[t-1, f] * (1 - x[t-1, f]);   // no ataque (decrece)

      theta[t, f] = inv_logit(alpha + beta * z[t, f]);
    }
  }

  // vectorizamos
  theta_vec = to_vector(theta);
  z_vec = to_vector(z);
}

model {
  // Previas sobre los parámetros, no sobre theta
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(beta_mu, beta_sigma);
  // a iota, delta y z1 les dejamos previas planas = uniformes en [0, 1].
  // No hace falta escribirlo.

  // Likelihood
  y ~ bernoulli(theta_vec);
}