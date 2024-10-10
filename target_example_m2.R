target <- 0 # log posterior
target <- target + dnorm(alpha, mean = alpha_mu, sd = alpha_sigma, log = T)
target <- target + dnorm(beta, mean = beta_mu, sd = beta_sigma, log = T)
target <- target + sum(dbinom(y, prob = theta, size = 1, log = T))