library(tidyverse)
library(markovchain)

normalize <- function(x) x / sum(x)

# Modelo recursivo con variable latente.
# La predictora observada (x = ataque) también sigue un modelo bernoulli
# recursivo.

# matriz de transición para x:
xmc <- new("markovchain", states = c("no attack", "attack"),
           transitionMatrix = xtmat <- matrix(c(0.8, 0.2,
                                                0.4, 0.6), 2, 2, byrow = T),
           name = "xmc")
(steady_x <- steadyStates(xmc))

# Para crear alphas y betas (intercepts y slopes) primero defino matrices de
# transición para z = 0 (tranquilas) y z = 1 (disturbadas)
# y = 1: en movimiento; y = 0: quietas.
ytmat0 <- new(
  "markovchain", states = c("quieta", "movimiento"),
  transitionMatrix = matrix(c(0.93, 0.07,
                              0.6, 0.4), 2, 2, byrow = T),
  name = "yz0"
)
steadyStates(ytmat0)

ytmat1 <- new(
  "markovchain", states = c("quieta", "movimiento"),
  transitionMatrix = matrix(c(0.55, 0.45,
                              0.2, 0.8), 2, 2, byrow = T),
  name = "yz1"
)
steadyStates(ytmat1)

# Definimos el modelo en base a los logits de p(y = 1: en movimiento)
alpha <- qlogis(ytmat0@transitionMatrix[, 2])
beta <- qlogis(ytmat1@transitionMatrix[, 2]) - qlogis(ytmat0@transitionMatrix[, 2])

curve(plogis(alpha[1] + beta[1] * x), ylim = c(0, 1), ylab = "P(y = 1)")
curve(plogis(alpha[2] + beta[2] * x), add = T, col = 2)

# check:
cbind(1 - plogis(alpha), plogis(alpha))
cbind(1 - plogis(alpha + beta), plogis(alpha + beta))

# dinámica de z
iota <- 0.7
delta <- 0.4

# Simulamos la dinámica de z para obtener estados iniciales, y de paso miramos
# la distribución estacionaria de y
set.seed(685)
nsim <- 3e4
use_start <- 1e4 + 1

xtmat <- xmc@transitionMatrix

zlong <- numeric(nsim)
xlong <- integer(nsim)
ylong <- integer(nsim)

zlong[1] <- runif(1)
xlong[1] <- rbinom(1, 1, steady_x[, 2])
ylong[1] <- rbinom(1, 1, steadyStates(ytmat0)[2])

for(i in 2:nsim) {
  # z: disturbio, en base al estado anterior y al ataque anterior
  zlong[i] <-
    zlong[i-1] +
    iota * (1 - zlong[i-1]) * xlong[i-1] -
    delta * zlong[i-1] * (1 - xlong[i-1])

  # y: comportamiento, en base al disturbio actual
  probs <- plogis(alpha + beta * zlong[i])
  ylong[i] <- rbinom(1, 1, probs[ylong[i-1] + 1])

  # x: ataque, sólo depende del ataque anterior
  xlong[i] <- rbinom(1, 1, xtmat[xlong[i-1] + 1, 2])
}

xs <- xlong[use_start:nsim]
zs <- zlong[use_start:nsim]
ys <- ylong[use_start:nsim]

table(xs) |> normalize()
table(ys) |> normalize()
hist(zs)

(tcont <- table(xs, ys))
tcont[, 2] / tcont[, 1]

dsim <- data.frame(y = ys, z = zs)
ggplot(dsim, aes(z, y)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ylim(0, 1)


# Simulamos datos

nt <- 12 # largo de cada seguimiento
nf <- 30 # nro de seguimientos

x <- y <- z <- matrix(NA, nt, nf)

z1 <- zs[sample(1:length(zs), nf, replace = F)]

# get steady probability of y[1] = 1 for each z0
steady_z1 <- sapply(1:nf, function(i) {
  mm <- cbind(
    1 - plogis(alpha + beta * z1[i]),
    plogis(alpha + beta * z1[i])
  )
  mmm <- new("markovchain", states = c("quieta", "movimiento"),
             transitionMatrix = mm,
             name = "transient")
  return(steadyStates(mmm)[2])
})

x[1, ] <- rbinom(nf, 1, steady_x[2])
y[1, ] <- rbinom(nf, 1, steady_z1)
z[1, ] <- z1

for(i in 2:nt) {
  # z (disturbio)
  z[i, ] <-
    z[i-1, ] +
    iota * (1 - z[i-1, ]) * x[i-1, ] -
    delta * z[i-1, ] * (1 - x[i-1, ])

  # y (comportamiento)
  for(f in 1:nf) {
    probs <- plogis(alpha + beta * z[i, f])
    y[i, f] <- rbinom(1, 1, probs[y[i-1, f] + 1])
  }

  # x (ataque)
  for(f in 1:nf) {
    x[i, f] <- rbinom(1, 1, xtmat[x[i-1, f] + 1, 2])
  }
}

table(x) |> normalize()
table(y) |> normalize()
tt <- table(x, y)
tt[, 2] / rowSums(tt)
# hist(z)

df <- data.frame(
  t = rep(1:nt, nf),
  grupo = rep(1:nf, each = nt),
  y = as.vector(y),
  x = as.vector(x),
  z = as.vector(z)
)

write.csv(df, "datos_ballenas.csv", row.names = F)

ggplot(df, aes(z, y)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ylim(0, 1)