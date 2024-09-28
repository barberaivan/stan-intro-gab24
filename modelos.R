# Paquetes ----------------------------------------------------------------

library(tidyverse); theme_set(theme_classic())
library(viridis)
library(rstan)
library(bbmle)
library(DHARMa)

# Funciones ---------------------------------------------------------------

# Para normalizar densidades o likelihoods (con fines gráficos).
# area es el tamaño del segmento o celda en la aproximación discreta.
# En 1D, es la distancia entre valores de la secuencia de parámetros.
# En 2D, es el área de cada celda.
# Devuelve un vector de densidad o likelihood, tal que
# sum(normalized_density * area) ~= 1
normalize_dens <- function(dens, area) {
  dens / sum(dens * area)
}

# resumir muestras de una distribución con media e intervalo de credibilidad,
# con probabilidad definida por ci
mean_ci <- function(x, ci = 0.95) {
  out <- 1 - ci
  qq <- quantile(x, probs = c(out / 2, ci + out / 2), method = 8) %>% unname
  return(c("mean" = mean(x), "lower" = qq[1], "upper" = qq[2]))
}

# logit y logit inverso de forma expresiva
logit <- function(x) qlogis(x)       # = log(x / (1 - x))
inv_logit <- function(x) plogis(x)   # = 1 / (1 + exp(-x))

# Tema elegante para ggplot
nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.35),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# Densidad y likelihood ---------------------------------------------------


# Datos y constantes ------------------------------------------------------

d <- read.csv("datos_ballenas.csv")
ncores <- min(parallel::detectCores() - 1, 4)

# Modelo 1: y ~ 1 ---------------------------------------------------------


# Likelihood -------------------------

# creamos una función que evalúe la likelihood para cualquier valor de p (theta).
# ll por log-likelihood
ll1 <- function(theta) {
  sum(dbinom(d$y, size = 1, prob = theta, log = T))
}

# versión vectorizada, para darle una secuencia de theta
ll1vec <- function(theta) {
  ll <- numeric(length(theta))
  for(i in 1:length(theta)) {
    ll[i] <- ll1(theta[i])
  }
  return(ll)
}

# graficamos para una secuencia fina de theta
theta_seq <- seq(0, 1, length.out = 10000)
ll_seq <- ll1vec(theta_seq)
plot(ll_seq ~ theta_seq, type = "l", ylab = "Log-Likelihood",
     xlab = expression(theta))
plot(exp(ll_seq) ~ theta_seq, type = "l", ylab = "Likelihood",
     xlab = expression(theta))

# Optimización a mano:
# el valor que tenga la mayor likelihood.
row_optim <- which.max(ll_seq)
(theta_hat1 <- theta_seq[row_optim])

# Optimización con bbmle::mle2.
# Creamos una función que devuelva la -log_likelihood, porque mle2 busca mínimos.
ll1_neg <- function(theta) -ll1(theta)

m1_mle <- mle2(ll1_neg,              # función de -log_likelihood
               list(theta = 0.5),    # valores iniciales
               method = "Brent",     # método de optimización, para 1D
               lower = 0, upper = 1) # límites del parámetro
summary(m1_mle)
coef(m1_mle)
bbmle::confint(m1_mle, level = 0.95)

# ploteamos
plot(exp(ll_seq) ~ theta_seq, type = "l", ylab = "Likelihood",
     xlab = expression(theta), xlim = c(0.2, 0.6))
# estimación a mano vs. la de mle2
abline(v = c(theta_hat1, coef(m1_mle)), col = c("gray60", "red"),
       lty = 1:2, lwd = c(5, 1))


# Posterior -------------------------

# Asumiendo una previa plana, la posterior = likelihood * previa
# será proporcional a esa likelihood. Ahora usaremos
# log(posterior) = log(likelihood) + log(previa)
# y la llamaremos "lp"

lp_un <- ll1vec(theta_seq) + dunif(theta_seq, log = T)
# "un" por unnomalized, o sea, una posterior que no integra 1.

par(mfrow = c(1, 2))
plot(ll_seq ~ theta_seq, type = "l", ylab = "Log-Likelihood",
     xlab = expression(theta))
plot(lp_un ~ theta_seq, type = "l", ylab = "Log-Posterior no normalizada",
     xlab = expression(theta))
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(exp(ll_seq) ~ theta_seq, type = "l", ylab = "Likelihood",
     xlab = expression(theta))
plot(exp(lp_un) ~ theta_seq, type = "l", ylab = "Posterior no normalizada",
     xlab = expression(theta))
par(mfrow = c(1, 1))

# Acá la densidad de una uniforme (probabilidad previa) es siempre 1,
# y como log(1) = 0, la likelihood y la posterior son iguales.


# Tomamos muestras de la posterior a mano, usando una aproximación discreta.
# Representamos la posterior con un set finito de npost valores, cada uno
# con probabilidad proporcional a la posterior no normalizada. Muestreamos
# ese set con reemplazo.
npost <- 10000
ids_post <- sample(1:length(theta_seq), size = npost, replace = T,
                   prob = exp(lp_un))
theta_draws <- theta_seq[ids_post]

# Ahora ploteamos la posterior normalizada, aproximada por la secuencia y
# la posterior aproximada por muestras

# normalizamos la secuencia
area <- diff(theta_seq[1:2])
dpost <- normalize_dens(exp(lp_un), area)

# Posterior analítica de la aproximación discreta
plot(dpost ~ theta_seq, type = "l", ylab = "Posterior normalizada",
     xlab = expression(theta), lwd = 5, col = "gray60")

# Densidad empírica de la posterior aproximada con muestras
lines(density(theta_draws, from = 0, to = 1, n = 2^10), col = "blue", lty = 2)

# Previa
lines(dunif(theta_seq) ~ theta_seq, col = "red")

# Usando muestras podemos mirar muchos detalles de la posterior, como
mean_ci(theta_draws, ci = 0.80) # media e intervalo del 80 %
mean_ci(theta_draws, ci = 0.95) # media e intervalo del 95 %

# Pr(theta > 0.3):
sum(theta_draws > 0.3) / npost

# Pr(0.3 < theta < 0.32)
sum(theta_draws > 0.3 & theta_draws < 0.32) / npost


# Muestremos la posterior con Stan ------------

# compilamos el modelo
smodel1 <- stan_model("modelo1.stan", verbose = F)

# preparamos datos
sdata1 <- list(n = nrow(d), y = d$y)

# muestreamos
m1 <- sampling(
  smodel1, data = sdata1, seed = 9663,
  cores = ncores, chains = 4, iter = 2000, warmup = 1000
)

# Chequeamos convergencia y calidad de las muestras
sm1 <- summary(m1)[[1]]
sm1["theta", c("n_eff", "Rhat")]
# Rhat debe ser < 1.01
# n_eff debe ser al menos > 400, pero idealmente > 1000

# Extraemos muestras y graficamos la posterior
theta_draws_stan <- as.matrix(m1)[, "theta", drop = T]

# Las mismas posteriores de antes
plot(dpost ~ theta_seq, type = "l", ylab = "Posterior normalizada",
     xlab = expression(theta), lwd = 5, col = "gray60")
lines(density(theta_draws, from = 0, to = 1, n = 2^10), col = "blue", lty = 2)
# Lo que muestreó Stan:
lines(density(theta_draws_stan, from = 0, to = 1, n = 2^10), col = "green", lty = 2)


# Comparamos estimación de Stan con MLE
c(coef(m1_mle), bbmle::confint(m1_mle, level = 0.95))
mean_ci(theta_draws_stan, ci = 0.95)


# Exploramos el efecto de la previa ---------------

# Aprovechando que en 1D podemos calcular la posterior casi analíticamente,
# exploramos el efecto de distintas previas.

# Usaremos una previa Beta, que está definida en [0, 1].
# La previa es una distribución de probabilidad bien conocida, por lo que
# estará normalizada (su área integra 1). Pero la likelihood no lo está,
# entonces generalmente no podemos comparar la forma de la previa, la likelihood
# y la posterior, a menos que normalicemos todo. Para eso usaremos la función
# normalize_dens.

# parámetros de la previa (para jugar, pero que ambos sean positivos)
a <- 1
b <- 6

# likelihood, prior, posterior
prior <- dbeta(theta_seq, a, b)  # normalizada
log_like <- ll1vec(theta_seq)    # no normalizada, porque no es una
                                 #  distribución
log_post <- log(prior) + log_like # no normalizada aún.

# Convertimos a escala de probabilidad (a partir de la log-probabilidad) y
# normalizamos
like <- normalize_dens(exp(log_like), area)
post <- normalize_dens(exp(log_post), area)

# Graficamos
ymax <- max(c(like, post, prior)) * 1.05
plot(prior ~ theta_seq, type = "l", col = "red", ylim = c(0, ymax),
     ylab = "Densidad o likelihood", xlab = expression(theta))
lines(like ~ theta_seq, col = "blue")
lines(post ~ theta_seq, col = "black")
text(0.9, ymax * 0.95, "Previa", col = "red")
text(0.9, ymax * 0.85, "Likelihood", col = "blue")
text(0.9, ymax * 0.75, "Posterior", col = "black")


# Modelo 2: y ~ z -------------------------------------------------------

# Fingimos haber observado la predictora latente, z, para tener un ejemplo
# de regresión logística.

# Función de likelihood
ll2 <- function(alpha, beta) {
  # calculamos theta en base a alpha, beta y z. Ahora theta no es un parámetro,
  # sino una cantidad derivada
  theta <- inv_logit(alpha + beta * d$z) # función logística, que restringe a
                                         # theta en [0, 1]
  ll <- sum(dbinom(d$y, size = 1, prob = theta, log = T))
  return(ll)
}

# versión vectorizada, para darle vectores de alpha y beta
ll2vec <- function(alpha, beta) {
  ll <- numeric(length(alpha))
  for(i in 1:length(alpha)) {
    ll[i] <- ll2(alpha[i], beta[i])
  }
  return(ll)
}

# Ahora el gráfico será 3D (o algo así), porque tenemos la likelihood
# en función de 2 parámetros.
# Para ello creamos una grilla de alpha y beta, pero lo ordenamos en una tabla,
# usando expand.grid
nside <- 250 # número de valores por parámetro. OJO: la grilla tendrá nside ^ 2 valores

# El rango de valores debería abarcar todo lo que sería razonable.
alpha_seq <- seq(-4, 4, length.out = nside)
beta_seq <- seq(-4, 4, length.out = nside)

partable <- expand.grid(alpha = alpha_seq,
                        beta = beta_seq)

partable$log_like <- ll2vec(partable$alpha, partable$beta)

# normalizamos, para luego poder comparar con la previa y la posterior.
area <- diff(alpha_seq[1:2]) * diff(beta_seq[1:2])
partable$like_norm <- normalize_dens(exp(partable$log_like), area = area)

ggplot(partable, aes(x = alpha, y = beta, fill = like_norm)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(fill = "Likelihood") +
  xlab(expression(alpha)) +
  ylab(expression(beta))

# Optimización a mano:

# buscamos la mayor likelihood en la grilla (partable, que es una grilla en
# forma de tabla)
row_optim <- which.max(partable$log_like)
(coef_hat1 <- partable[row_optim, c("alpha", "beta")])

# Optimización con mle2
ll2_neg <- function(alpha, beta) -ll2(alpha, beta)

m2_mle <- mle2(ll2_neg,                   # función de -log_likelihood
               start = list(alpha = 0, beta = 0), # valores iniciales
               lower = list(alpha = -Inf, beta = -Inf),
               upper = list(alpha = Inf, beta = Inf),
               method = "L-BFGS-B") # Este método permite poner límites en más de 1D.
# En este caso es absurdo poner límites, pero sirve de ejemplo.

summary(m2_mle) # HAY ESTRELLITAS
coef(m2_mle)
bbmle::confint(m2_mle, level = 0.95)


# Metemos ambas estimaciones en un data.frame para ggplot
coef_hat_table <- data.frame(
  alpha = c(coef_hat1[1, 1], coef(m2_mle)[1]),
  beta = c(coef_hat1[1, 2], coef(m2_mle)[2]),
  metodo = c("grilla", "mle2")
)

# ploteamos nuevamente, pero haciendo zoom sobre la zona no-nula de la
# likelihood e incluyendo las estimaciones puntuales
ggplot(partable, aes(x = alpha, y = beta, fill = like_norm)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(fill = "Likelihood") +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  # agregamos las estimaciones puntuales
  geom_point(data = coef_hat_table, aes(alpha, beta, shape = metodo),
             inherit.aes = F, size = 4) +
  scale_shape_manual(values = c(1, 4)) +
  xlim(-3, -1) + # esto hace zoom
  ylim(1.5, 3.8)


# Likelihood, previa y posterior -------------

# Graficaremos las 3 variables en este formato, para explorar el efecto de
# la previa.
# Como alpha y beta no tienen restricciones, podemos usar previas normales.

# Parámetros de las previas
alpha_mu <- 0
alpha_sigma <- 0.25

beta_mu <- 2
beta_sigma <- 0.25

# Previa
partable$log_prior <-
  dnorm(partable$alpha, mean = alpha_mu, sd = alpha_sigma, log = T) +
  dnorm(partable$beta, mean = beta_mu, sd = beta_sigma, log = T)
partable$prior <- exp(partable$log_prior)

# Posterior
partable$log_post <- partable$log_like + partable$log_prior
partable$post <- normalize_dens(exp(partable$log_post), area)

# Ordenamos para plotear
nn <- which(names(partable) %in% c("like_norm", "prior", "post"))
ptlong <- pivot_longer(partable,
                       cols = all_of(nn), names_to = "var",
                       values_to = "densidad")
ptlong$var <- factor(ptlong$var, levels = c("like_norm", "prior", "post"),
                     labels = c("Likelihood", "Previa", "Posterior"))

ggplot(ptlong, aes(x = alpha, y = beta, fill = densidad)) +
  geom_tile(na.rm = T) +
  scale_fill_viridis(na.value = "transparent") +
  labs(fill = "Densidad") +
  facet_wrap(vars(var), nrow = 2, axes = "all", axis.labels = "margins") +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  theme(strip.background = element_rect(color = "white"),
        legend.position.inside = c(0.75, 0.25)) +
  xlim(-3, 1.75) +
  ylim(1, 3.8)


# Muestreamos la posterior a mano -------------

# Volvemos a definir previas, así nos independizamos de las elecciones que
# hicimos para el gráfico. También hay que recalcular la posterior

alpha_mu <- 0
alpha_sigma <- 5
beta_mu <- 0
beta_sigma <- 2

# Previa
partable$log_prior <-
  dnorm(partable$alpha, mean = alpha_mu, sd = alpha_sigma, log = T) +
  dnorm(partable$beta, mean = beta_mu, sd = beta_sigma, log = T)
partable$prior <- exp(partable$log_prior)

# Posterior
partable$log_post <- partable$log_like + partable$log_prior
partable$post <- normalize_dens(exp(partable$log_post), area)

# Igual que antes, sorteamos filas de partable según la probabilidad posterior.
npost <- 10000
ids_post <- sample(1:nrow(partable), size = npost, replace = T,
                   prob = partable$post)
m2_draws_hand <- as.matrix(partable[ids_post, c("alpha", "beta")])

# Muestreamos la posterior con Stan -------------

# compilamos el modelo
smodel2 <- stan_model("modelo2.stan", verbose = F)

# preparamos datos
sdata2 <- list(
  n = nrow(d), y = d$y, z = d$z,
  alpha_mu = alpha_mu,
  alpha_sigma = alpha_sigma,
  beta_mu = beta_mu,
  beta_sigma = beta_sigma
)

# muestreamos
m2 <- sampling(
  smodel2, data = sdata2, seed = 9663,
  cores = ncores, chains = 4, iter = 2000, warmup = 1000
)

# Chequeamos convergencia y calidad de las muestras
sm2 <- summary(m2)[[1]]
sm2[c("alpha", "beta"), c("n_eff", "Rhat")]

# Extraemos muestras y graficamos la posterior
m2_draws <- as.matrix(m2, pars = c("alpha", "beta"))

# Comparamos muestreo a mano con el de Stan

# Distribuciones marginales
par(mfrow = c(1, 2))
plot(density(m2_draws[, "alpha"]), main = NA, xlab = expression(alpha))
lines(density(m2_draws_hand[, "alpha"]), col = "blue")

plot(density(m2_draws[, "beta"]), main = NA, xlab = expression(beta))
lines(density(m2_draws_hand[, "beta"]), col = "blue")
par(mfrow = c(1, 1))

# Distribución conjunta
xr <- range(rbind(m2_draws, m2_draws_hand)[, "beta"])
xr[1] <- xr[1] - abs(diff(xr)) * 0.05
xr[2] <- xr[2] + abs(diff(xr)) * 0.05

yr <- range(rbind(m2_draws, m2_draws_hand)[, "alpha"])
yr[1] <- yr[1] - abs(diff(yr)) * 0.05
yr[2] <- yr[2] + abs(diff(yr)) * 0.05

par(mfrow = c(1, 2))
plot(m2_draws[, "alpha"] ~ m2_draws[, "beta"], main = "Stan",
     ylab = expression(alpha), xlab = expression(beta),
     col = rgb(0, 0, 0, 0.1), pch = 19, ylim = yr, xlim = xr)

plot(m2_draws_hand[, "alpha"] ~ m2_draws_hand[, "beta"], main = "A mano",
     ylab = expression(alpha), xlab = expression(beta),
     col = rgb(0, 0, 0, 0.1), pch = 19, ylim = yr, xlim = xr)
par(mfrow = c(1, 1))


# Evaluación de modelos con DHARMa ----------------

# Los modelos estadísticos permiten simular datos. DHARMa calcula residuos de
# probabilidad acumulada aleatorizada, que bajo un buen ajuste tienen
# distribución uniforme. Para calcularlos se requieren los datos observados
# (respuesta) y una matriz de datos simulados.

nsim <- 1000 # nro de sets de datos a simular

# Modelo estimado por maximum likelihood.
coef_mle <- coef(m2_mle)
ysim_mle <- matrix(NA, nrow(d), nsim)
theta_fitted <- inv_logit(coef_mle["alpha"] + coef_mle["beta"] * d$z)
for(i in 1:nsim) {
  ysim_mle[, i] <- rbinom(n = nrow(d), size = 1, prob = theta_fitted)
}

# Modelo estimado con enfoque Bayesiano, con Stan
ysim_stan <- matrix(NA, nrow(d), nsim)
# para considerar incertidumbre, cada vez que simulamos datos tomaremos una
# muestra de la posterior
npost_stan <- nrow(m2_draws)
# sorteamos 1000 muestras previamente
ids_res <- sample(1:npost_stan, size = nsim, replace = F)
for(i in 1:nsim) {
  row <- ids_res[i]
  theta_fitted_i <- inv_logit(
    m2_draws[row, "alpha"] +
    m2_draws[row, "beta"] * d$z
  )
  ysim_stan[, i] <- rbinom(n = nrow(d), size = 1, prob = theta_fitted_i)
}

# Calculamos residuos
res2_mle <- createDHARMa(observedResponse = d$y,
                         simulatedResponse = ysim_mle)
res2_stan <- createDHARMa(observedResponse = d$y,
                          simulatedResponse = ysim_stan)

plot(res2_mle)
plot(res2_stan)

# Residuos en función de la predictora (z)
plotResiduals(res2_mle, form = d$z, rank = F)
plotResiduals(res2_stan, form = d$z, rank = F)


# Graficamos la función ajustada (predicciones) ---------

# Arriba calculamos la probabilidad estimada para cada dato (theta_fitted),
# usando el valor observado de la predictora. Para visualizar la función
# ajustada por el modelo es más conveniente generar una secuencia de valores de
# z (predictora) y evaluar theta en cada valor.

nseq <- 200
zseq <- seq(0, 1, length.out = nseq)

# La predicción basada en la estimación puntual es la siguiente:
theta_pred_mle <- inv_logit(coef_mle["alpha"] + coef_mle["beta"] * zseq)
plot(theta_pred_mle ~ zseq, type = "l", ylab = expression(theta), xlab = "z",
     ylim = c(0, 1))

# Pero eso no basta, necesitamos un intervalo de confianza. Una forma de obtener
# un intervalo aproximado es aproximar la incertidumbre con una distribución
# normal, esa misma con la que obtenemos el "Std. Error" de los coeficientes
# cuando hacemos
summary(m2_mle)
# Pero la aproximación normal es válida sólo en la escala lineal, o sea, en la
# escala logit (el predictor lineal en la jerga de GLM):
theta_logit <- coef_mle["alpha"] + coef_mle["beta"] * zseq
# Calculando el standard error de la aproximación lineal podemos obtener el IC, y
# luego transformar a la escala de la probabilidad (inv_logit).

V <- vcov(m2_mle)
X <- cbind(rep(1, nseq), zseq) # matriz de diseño
theta_logit_se <- sqrt(diag(X %*% V %*% t(X))) # standard error de la predicción.
# Esto es lo que hace la función predict(..., se.fit = T) cuando le pasamos un
# GLM.

# Ponemos todo en un data.frame
pred_mle <- data.frame(
  z = zseq,
  theta = theta_pred_mle,
  theta_lower = inv_logit(theta_logit - qnorm(0.975) * theta_logit_se),
  theta_upper = inv_logit(theta_logit + qnorm(0.975) * theta_logit_se)
)

ggplot(pred_mle, aes(z, theta, ymin = theta_lower, ymax = theta_upper)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  ylim(0, 1) +
  ylab(expression(theta))

# Con el enfoque bayesiano hay que hacer un poco más de lío, pero es
# conceptualmente más simple: calcular lo deseado con cada muestra de la
# posterior, y al final de todo, resumir el resultado para plotear.

# Como tenemos 4000 muestras, generaremos una matriz de predicciones con 4000
# columnas
theta_pred <- matrix(NA, nseq, npost_stan)
# acá hacemos la misma cuenta para theta que hicimos arriba para simular datos,
# pero usando zseq en vez de d$z
for(i in 1:npost_stan) {
  theta_pred[, i] <- inv_logit(m2_draws[i, "alpha"] +
                               m2_draws[i, "beta"] * zseq)
}
# cada fila de theta pred contiene muestras de la distribución posterior de theta
# predicho para cada valor en zseq.
plot(density(theta_pred[1, ]))
plot(density(theta_pred[50, ]))

# Cada columna tiene una función, correspondiente a los valores de alpha y beta
# de cada muestra. Graficamos algunas
plot(theta_pred[, 1] ~ zseq, ylim = c(0, 1), ylab = expression(theta),
     xlab = "z", col = rgb(0, 0, 0, 0.2), type = "l")
for(i in sample(1:npost_stan, 200)) {
  lines(theta_pred[, i] ~ zseq, col = rgb(0, 0, 0, 0.2))
}

# Y ahora resumimos la posterior de cada fila para hacer un gráfico de media
# e intervalo de credibilidad.
pred_stan <- as.data.frame(cbind(
 z = zseq,
 apply(theta_pred, 1, mean_ci) |> t() # esto resume con media e ic
))

ggplot(pred_stan, aes(z, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  ylim(0, 1) +
  ylab(expression(theta))

# Comparamos ambas predicciones
colnames(pred_stan) <- colnames(pred_mle)
pred_stan$approach <- "Bayesiano"
pred_mle$approach <- "Frecuentista"

pred_bayentista <- rbind(pred_stan, pred_mle)
ggplot(pred_bayentista,
       aes(z, theta, ymin = theta_lower, ymax = theta_upper)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(vars(approach), nrow = 1, axes = "all", axis.labels = "margins") +
  ylim(0, 1) +
  ylab(expression(theta)) +
  nice_theme()



# Modelo 3: y ~ z, z ~ x, con z latente -----------------------------------

# Ahora los parámetros son

# alpha
# beta
# iota
# delta
# z1 (vector de largo nf)

# En tantas dimensiones no podremos hacer gráficos como antes, así que
# directamente ajustaremos los modelos con mle2 y con stan.

## Ordenamos los datos en forma matricial
nt <- max(d$t)
nf <- max(d$grupo) # f de follow, seguimiento.

y <- matrix(d$y, nt, nf)
x <- matrix(d$x, nt, nf)
# z es latente! no la vimos.

# Función de likelihood
ll3 <- function(params) {
  # Extraemos los parámetros del vector params
  alpha <- params[1]
  beta <- params[2]
  iota <- params[3]
  delta <- params[4]
  z1 <- params[5:(nf+4)]

  # matriz de predictora latente
  z <- matrix(NA, nt, nf)
  z[1, ] <- z1           # primer estado se estima

  theta <- matrix(NA, nt, nf) # cantidad derivada
  theta[1, ] <- inv_logit(alpha + beta * z[1, ])

  # Loop recursivo para calcular theta y z
  for(t in 2:nt) {
    # actualizamos z en base al ataque en el intervalo anterior
    z[t, ] <-
      z[t-1, ] +
      iota * (1 - z[t-1, ]) * x[t-1, ] -  # ataque (aumenta)
      delta * z[t-1, ] * (1 - x[t-1, ])   # no ataque (decrece)

    theta[t, ] <- inv_logit(alpha + beta * z[t, ])
  }

  # vectorizamos theta, para calcular likelihood.
  ll <- sum(dbinom(d$y, size = 1, prob = as.vector(theta), log = T))
  return(ll)
}

# La negamos para mle2
ll3_neg <- function(x) -ll3(x)

# nombres
parnames1 <- c("alpha", "beta", "iota", "delta")
parnames2 <- paste("z1", 1:nf, sep = "_")
parnames <- c(parnames1, parnames2)

# valores iniciales
start1 <- c(0, 1, 0.5, 0.5)
start2 <- rep(0.5, nf) # para z1
start <- c(start1, start2)
names(start) <- parnames

# bounds
lower1 <- c(-Inf, -Inf, 0, 0)
upper1 <- c(Inf, Inf, 1, 1)
lower2 <- rep(0, nf)
upper2 <- rep(1, nf)

lower <- c(lower1, lower2)
upper <- c(upper1, upper2)
names(lower) <- names(upper) <- parnames

# Estimamos
parnames(ll3_neg) <- parnames
m3_mle <- mle2(ll3_neg, start = start, parnames = parnames, vecpar = T,
               lower = lower, upper = upper, method = "L-BFGS-B")
summary(m3_mle) # Hay problemas. Mejor ser bayesianos.


## Muestreamos la posterior con Stan

smodel3 <- stan_model("modelo3.stan", verbose = F)

# preparamos datos
sdata3 <- list(
  nt = nt, nf = nf, y = d$y, x = x,

  # previas
  alpha_mu = alpha_mu,
  alpha_sigma = alpha_sigma,
  beta_mu = beta_mu,
  beta_sigma = beta_sigma
)

# muestreamos
m3 <- sampling(
  smodel3, data = sdata3, seed = 9663,
  cores = ncores, chains = 4, iter = 2000, warmup = 1000,
  pars = c("alpha", "beta", "iota", "delta", "z1",
           "z_vec", "theta_vec")
)

# Chequeamos convergencia y calidad de las muestras
sm3 <- summary(m3)[[1]]
min(sm3[, "n_eff"])
max(sm3[, "Rhat"])
# Biutiful

# Miramos rápidamente las marginales, y comparamos con las previas
m3_draws <- as.matrix(m3)

plot(density(m3_draws[, "alpha"]), main = NA, xlab = expression(alpha))
curve(dnorm(x, alpha_mu, alpha_sigma), add = T, col = 2)

plot(density(m3_draws[, "beta"]), main = NA, xlab = expression(beta))
curve(dnorm(x, beta_mu, beta_sigma), add = T, col = 2)

plot(density(m3_draws[, "iota"], from = 0, to = 1),
     main = NA, xlab = expression(iota))
curve(dunif(x, 0, 1), add = T, col = 2)

plot(density(m3_draws[, "delta"], from = 0, to = 1),
     main = NA, xlab = expression(delta))
curve(dunif(x, 0, 1), add = T, col = 2)

# Los z1
z1_draws <- m3_draws[, grep("z1", colnames(m3_draws))]
colorcines <- viridis(nf, option = "H")

plot(density(z1_draws[, 1], from = 0, to = 1, adjust = 1.5), main = NA, xlab = "z1",
     ylim = c(0, 4.5), col = colorcines[1])
for(i in 2:nf) {
  lines(density(z1_draws[, i], from = 0, to = 1, adjust = 1.5),
        col = colorcines[i])
}
curve(dunif(x, 0, 1), add = T, col = "black", lwd = 1.5, lty = 2)


# ¿Y los z en t = nt?
zvec_draws <- m3_draws[, grep("z_vec", colnames(m3_draws))]
colnames(zvec_draws)
cols_nt <- d$t == nt
znt <- zvec_draws[, cols_nt] # sólo los z para t = nt

plot(density(znt[, 1], from = 0, to = 1, adjust = 1.5), main = NA, xlab = "z1",
     ylim = c(0, 35), col = colorcines[1])
for(i in 2:nf) {
  lines(density(znt[, i], from = 0, to = 1, adjust = 1.5),
        col = colorcines[i])
}
curve(dunif(x, 0, 1), add = T, col = "black", lwd = 1.5, lty = 2)


## Volvemos a intentar la likelihood

# Pero ahora estimamos todos los parámetros en la escala no restringida,
# es decir, los transformamos para que puedan estar en (-Inf, Inf).
# A la vez, para evitar que se vayan al infinito, usamos una likelihood
# penalizada.

# La estimación de máxima verosimilitud restringida no minimiza la
# -log_likelihood sino que define como objetivo
# -log_likelihood + penalidad
# Y la penalidad se diseña de tal forma que no permita que el mínimo de
# ese término esté en +-Inf.
# Por ejemplo, el método de regularización y selección de variables LASSO
# (Least Angular Shrinkage and Selection Operator) usa la siguiente penalidad:
# -log_likelihood + lambda * sum(beta ^ 2)
# Donde beta es el vector de coeficientes de regresión, y lambda > 0
# controla el grado de penalización. De esta manera, para valores razonablemente
# altos de lambda, la penalidad aumenta con el valor absoluto de los betas,
# y no podemos tener un mínimo (óptimo) en beta = +-Inf.

# la -log_posterior tiene una forma similar: -log_likelihood + (-log_previa).
# Si minimizamos la -log_posterior o maximizamos la log_posterior obtenemos
# el modo de la posterior, conocido como MAP, por Maximum a Posteriori.
# En este caso, sería como minimizar la -log_likelihood penalizada,
# donde la penalidad es -log_previa.

# Entonces podemos pensar en estimar una likelihood penalizada o estimar el MAP
# de una posterior.

# Los parámetros iota, delta y z1 deben estar en [0, 1]. Usaremos sus logits
# para que vivan en la escala (-Inf, Inf). Y como queremos asignarles una previa
# uniforme en su escala original, les asignaremos una previa logística en la
# escala logit. Sí, magia, funciona:

samples_unif <- runif(1e5, 0, 1)
samples_logit <- logit(samples_unif)
plot(density(samples_logit), lwd = 1.5)
curve(dlogis(x), add = T, col = "red", lty = 2, lwd = 1.5)

# Esto evitará que el algoritmo busque el óptimo en +-Inf.
# A la vez, aproximar la incertidumbre con una normal multivariada será más
# razonable si no estamos contra el borde del espacio de parámetros.

# Función de likelihood penalizada (similar a una posterior)
ll3_pen <- function(params) {
  # Extraemos los parámetros del vector params
  alpha <- params[1]
  beta <- params[2]       # alpha y beta sin penalizar
  iota_logit <- params[3]
  delta_logit <- params[4]
  z1_logit <- params[5:(nf+4)]

  # pasamos a la escala original para calcular la likelihood
  iota <- inv_logit(iota_logit)
  delta <- inv_logit(delta_logit)
  z1 <- inv_logit(z1_logit)

  # matriz de predictora latente
  z <- matrix(NA, nt, nf)
  z[1, ] <- z1           # primer estado se estima

  theta <- matrix(NA, nt, nf) # cantidad derivada
  theta[1, ] <- inv_logit(alpha + beta * z[1, ])

  # Loop recursivo para calcular theta y z
  for(t in 2:nt) {
    # actualizamos z en base al ataque en el intervalo anterior
    z[t, ] <-
      z[t-1, ] +
      iota * (1 - z[t-1, ]) * x[t-1, ] -  # ataque (aumenta)
      delta * z[t-1, ] * (1 - x[t-1, ])   # no ataque (decrece)

    theta[t, ] <- inv_logit(alpha + beta * z[t, ])
  }

  # vectorizamos theta, para calcular likelihood.
  ll <- sum(dbinom(d$y, size = 1, prob = as.vector(theta), log = T))

  # Ahora calculamos la log_previa, que sería el término de penalización
  # pero sólo asignamos previa a los parámetros complicados:
  # iota_logit, delta_logit y z1_logit
  log_prior <-
    dlogis(iota_logit, log = T) + dlogis(delta_logit, log = T) +
    sum(dlogis(z1_logit, log = T))

  return(ll + log_prior)
}

# La negamos para mle2
ll3_pen_neg <- function(x) -ll3_pen(x)

# valores iniciales
start1 <- c(0, 1, 0, 0) # iota_logit y delta_logit empiezan en 0, que implica
                        # iota = delta = 0.5
start2 <- rep(0, nf)    # lo mismo para para z1_logit
start <- c(start1, start2)
names(start) <- parnames

# Estimamos
parnames(ll3_pen_neg) <- parnames

m3_mle_pen <- mle2(ll3_pen_neg, start = start, parnames = parnames, vecpar = T,
                   method = "BFGS")
# No le damos lower y upper, porque ahora ningún parámetro está restringido.
summary(m3_mle_pen) # Ahora tenemos el std.error para todos los parámetros.

cc <- coef(m3_mle_pen)
inv_logit(cc["iota"])
inv_logit(cc["delta"])
inv_logit(cc[grep("z1", names(cc))]) |> plot()
