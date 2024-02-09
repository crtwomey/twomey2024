#
# rd.R
#
# Reference implementation of the Blahut-Arimoto algorithm
# for rate-distortion quantization.
#
# Blahut-Arimoto with an Expectation-Maximization-like step
# for estimating the optimal quantization centers. Achieves
# a local optimum for Bregman divergences (e.g. squared
# Euclidean distance and KL-divergence).
#
# All code is open source under the GPL v3 (see LICENSE)
#

# all-pairs distances
squared_euclidean_distance <- function(x, y) {
  s <- matrix(0, nrow(x), nrow(y))
  for (i in seq_len(ncol(x))) {
    s <- s + outer(x[, i], y[, i], "-")^2
  }
  return(s)
}

# all-pairs KL-divergences
kl_divergence <- function(px, qx) {
  s <- matrix(0, nrow(px), nrow(qx))
  for (i in seq_len(ncol(px))) {
    s <- s + px[, i] * outer(log(px[, i]), log(qx[, i]), "-")
  }
  return(s)
}

fqxhat_x <- function(beta, pxhat, x, xhat, d) {
  pxhat * exp(-beta * t(d(x, xhat)))
}

fpxhat_x <- function(beta, pxhat, x, xhat, d) {
  qxhat_x <- fqxhat_x(beta, pxhat, x, xhat, d)
  qxhat_x <- qxhat_x + 1E-16
  zxhat   <- colSums(qxhat_x)
  sweep(qxhat_x, 2, zxhat, "/")
}

# Sample initial centers from x
init_rd <- function(x,
  nxhat = 5,
  beta  = 50,
  px    = matrix(1 / nrow(x), nrow(x), 1),
  fd    = squared_euclidean_distance,
  xhat  = c()
) {
  if (length(xhat) == 0) {
    s    <- sample(seq_len(nrow(x)), nxhat)
    xhat <- as.matrix(x[s, ], nxhat, ncol(x))
  } else {
    nxhat <- nrow(xhat)
  }
  nx      <- length(px)
  pxhat   <- rep(1 / nxhat, nxhat)
  pxhat_x <- matrix(1 / nxhat, nxhat, nx)
  list(
    nxhat   = nxhat,
    beta    = beta,
    d       = fd,
    pxhat_x = pxhat_x,
    pxhat   = pxhat,
    xhat    = xhat,
    x       = x,
    px      = px
  )
}

# uses different sampling strategy for initialization
init_dirichlet_rd <- function(x,
  nxhat = 5,
  beta  = 50,
  px    = matrix(1 / nrow(x), nrow(x), 1),
  fd    = squared_euclidean_distance,
  xhat  = c()
) {
  if (length(xhat) == 0) {
    s    <- sample(seq_len(nrow(x)), nxhat)
    xhat <- as.matrix(x[s, ], nxhat, ncol(x))
  } else {
    nxhat <- nrow(xhat)
  }
  nx      <- length(px)
  pxhat   <- rep(1 / nxhat, nxhat)
  pxhat_x <- t(gtools::rdirichlet(nx, rep(1, nxhat)))
  list(
    nxhat   = nxhat,
    beta    = beta,
    d       = fd,
    pxhat_x = pxhat_x,
    pxhat   = pxhat,
    xhat    = xhat,
    x       = x,
    px      = px
  )
}

# Blahut-Arimoto iterative update
update_rd <- function(rd,
  update_xhat  = TRUE,
  update_pxhat = TRUE
) {
  pxhat_x <- fpxhat_x(rd$beta, rd$pxhat, rd$x, rd$xhat, rd$d)
  pxhat   <- rd$pxhat
  if (update_pxhat) {
    pxhat <- as.vector(pxhat_x %*% rd$px)
  }
  px_xhat <- sweep(pxhat_x, 2, rd$px, "*")
  px_xhat <- sweep(px_xhat, 1, rowSums(px_xhat), "/")
  xhat    <- rd$xhat

  if (update_xhat) {
    xhat <- px_xhat %*% rd$x
  }

  list(
    nxhat   = rd$nxhat,
    beta    = rd$beta,
    d       = rd$d,
    pxhat_x = pxhat_x,
    pxhat   = pxhat,
    xhat    = xhat,
    x       = rd$x,
    px      = rd$px
  )
}

# computes the average distortion for a given rd solution
average_distortion <- function(rd) {
  dx_xhat <- t(rd$d(rd$x, rd$xhat))
  pxhatx <- rd$pxhat_x * rd$px
  sum(pxhatx * dx_xhat)
}

# computes the mutual information for a given rd solution
mutual_information <- function(rd) {
  pxhatx <- sweep(rd$pxhat_x, 2, rd$px, "*")
  sum(pxhatx * log2(pxhatx / outer(rowSums(pxhatx), colSums(pxhatx), "*")))
}

# RMSE between centroids for two vocabularies
delta_rd <- function(rd1, rd2) {
  sum(rowSums((rd1$xhat - rd2$xhat)^2))
}

# redefine based on best mapping between two vocabularies
best_mapped_delta_rd <- function(rd1, rd2) {
  xhat <- rd1$xhat
  yhat <- rd2$xhat
  n    <- nrow(xhat)
  m    <- nrow(yhat)
  if (n != m) {
    return(Inf)
  }
  cost <- squared_euclidean_distance(xhat, yhat)
  map  <- lpSolve::lp.assign(cost, direction = "min")
  return((1 / n) * map$objval)
}
