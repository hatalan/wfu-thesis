#!/usr/bin/env Rscript
# Data generating processes for the VAR lag-misspecification sweeps (var-sweep-02 on).
#
# Four bivariate DGPs, all sharing the baseline innovation covariance Sigma and mean mu:
#
#   baseline        VAR(2) from var-2-notebook.qmd (used in var-sweep-01).
#   long_decay      VAR(6) with lag-i coefficients proportional to A1 / i (Minnesota-style
#                   1/i decay), scaled so its population R2_trace equals the baseline's.
#   sparse          Only lags 1 and 4 are nonzero: lag 1 has the baseline A1 pattern and
#                   lag 4 a diagonal own-lag "seasonal" term, scaled to the baseline R2_trace.
#   near_unit_root  Baseline VAR(2) with A_i -> kappa^i A_i. This multiplies every companion
#                   eigenvalue by kappa, giving spectral radius 0.97 with the same dynamics
#                   shape. Its R2 is high by construction (persistence raises predictability).
#
# Usage
#   source("notebooks/var/var-dgp.R")      # from the repo root (or "var-dgp.R" in notebooks/var)
#   dgps <- make_dgps()                    # list of DGPs with A, Sigma, mu, c, summary stats
#   Y    <- dgp_simulate(dgps$sparse, T_obs = 200, seed = 1)
#
#   Rscript notebooks/var/var-dgp.R        # print a summary table of the DGPs

# ---- Population quantities --------------------------------------------------

# Companion matrix of a VAR(p); A is a list of N x N lag matrices.
dgp_companion <- function(A) {
  p <- length(A)
  N <- nrow(A[[1]])
  F <- matrix(0, N * p, N * p)
  F[1:N, ] <- do.call(cbind, A)
  if (p > 1) F[(N + 1):(N * p), 1:(N * (p - 1))] <- diag(N * (p - 1))
  F
}

# Largest companion eigenvalue modulus; < 1 means stationary.
dgp_spectral_radius <- function(A) {
  max(Mod(eigen(dgp_companion(A), only.values = TRUE)$values))
}

# Unconditional covariance Gamma0 of y_t: vec(G) = (I - F (x) F)^{-1} vec(Q).
dgp_gamma0 <- function(A, Sigma) {
  F  <- dgp_companion(A)
  Np <- nrow(F)
  N  <- nrow(Sigma)
  Q  <- matrix(0, Np, Np)
  Q[1:N, 1:N] <- Sigma
  G <- matrix(solve(diag(Np^2) - kronecker(F, F), as.vector(Q)), Np, Np)
  G[1:N, 1:N]
}

# Population R2 variants (same definitions as the Stan generated quantities).
dgp_r2 <- function(A, Sigma) {
  if (dgp_spectral_radius(A) >= 1) stop("VAR is not stationary")
  G  <- dgp_gamma0(A, Sigma)
  eq <- 1 - diag(Sigma) / diag(G)
  list(R2_eq = eq, R2_mean = mean(eq), R2_trace = 1 - sum(diag(Sigma)) / sum(diag(G)))
}

# ---- Constructors ------------------------------------------------------------

# Scale every lag matrix by a common s so that R2_trace hits `target`.
dgp_scale_to_r2 <- function(shape, Sigma, target) {
  scaled <- function(s) lapply(shape, `*`, s)
  f <- function(s) {
    A <- scaled(s)
    if (dgp_spectral_radius(A) >= 1) return(1 - target)  # beyond the boundary R2 -> 1
    dgp_r2(A, Sigma)$R2_trace - target
  }
  hi <- 1
  while (dgp_spectral_radius(scaled(hi)) < 1) hi <- 2 * hi
  s <- uniroot(f, c(1e-8, hi), tol = 1e-12)$root
  structure(scaled(s), scale = s)
}

# Multiply every companion eigenvalue by kappa via A_i -> kappa^i A_i,
# with kappa chosen to give spectral radius `target_radius`.
dgp_scale_radius <- function(A, target_radius) {
  kappa <- target_radius / dgp_spectral_radius(A)
  structure(Map(function(Ai, i) kappa^i * Ai, A, seq_along(A)), scale = kappa)
}

# ---- DGP definitions -----------------------------------------------------------

dgp_baseline_spec <- list(
  A1    = matrix(c(0.45, -0.15, -0.35, 0.30), nrow = 2),   # filled by column
  A2    = matrix(c(0.05, -0.05, -0.05, 0.03), nrow = 2),
  c_vec = c(2.5, 1.5),
  Sigma = matrix(c(1.0, 0.1, 0.1, 1.0), nrow = 2)
)

make_dgps <- function(spec = dgp_baseline_spec,
                      long_lags = 6, sparse_lag = 4, sparse_own = 0.5,
                      near_unit_radius = 0.97) {
  A1    <- spec$A1
  A2    <- spec$A2
  Sigma <- spec$Sigma
  N     <- nrow(A1)
  Z     <- matrix(0, N, N)

  A_base <- list(A1, A2)
  mu     <- solve(diag(N) - A1 - A2, spec$c_vec)   # baseline mean, shared by all DGPs
  target <- dgp_r2(A_base, Sigma)$R2_trace

  A_list <- list(
    baseline       = A_base,
    long_decay     = dgp_scale_to_r2(lapply(seq_len(long_lags), function(i) A1 / i),
                                     Sigma, target),
    sparse         = dgp_scale_to_r2(c(list(A1), rep(list(Z), sparse_lag - 2),
                                       list(sparse_own * diag(N))),
                                     Sigma, target),
    near_unit_root = dgp_scale_radius(A_base, near_unit_radius)
  )

  descriptions <- c(
    baseline       = "VAR(2) from var-2-notebook.qmd",
    long_decay     = sprintf("VAR(%d), A_i proportional to A1 / i, R2 matched to baseline", long_lags),
    sparse         = sprintf("Lags 1 and %d only, R2 matched to baseline", sparse_lag),
    near_unit_root = sprintf("Baseline with eigenvalues scaled to radius %.2f", near_unit_radius)
  )

  Map(function(A, name) {
    r2 <- dgp_r2(A, Sigma)
    list(
      name            = name,
      description     = unname(descriptions[name]),
      A               = unname(lapply(A, unname)),
      p               = length(A),
      Sigma           = Sigma,
      mu              = mu,
      c_vec           = as.vector((diag(N) - Reduce(`+`, A)) %*% mu),
      scale           = if (is.null(attr(A, "scale"))) 1 else attr(A, "scale"),
      spectral_radius = dgp_spectral_radius(A),
      R2_eq           = r2$R2_eq,
      R2_mean         = r2$R2_mean,
      R2_trace        = r2$R2_trace
    )
  }, A_list, names(A_list))
}

# One row per DGP: true lag order, persistence and population R2.
dgp_summary <- function(dgps) {
  do.call(rbind, lapply(dgps, function(d) {
    data.frame(dgp = d$name, p = d$p, scale = d$scale,
               spectral_radius = d$spectral_radius,
               R2_eq1 = d$R2_eq[1], R2_eq2 = d$R2_eq[2],
               R2_mean = d$R2_mean, R2_trace = d$R2_trace,
               row.names = NULL)
  }))
}

# ---- Block DGP (var-sweep-03) -------------------------------------------------

# An 8-variable VAR(2) made of four independent 2-variable blocks. Block k is the
# baseline VAR(2) with both lag matrices scaled so its population R2_trace equals
# block_r2[k]; Sigma is block-diagonal with the baseline Sigma in each block. Because
# there are no links between blocks, the first N = 2, 4, 6 or 8 series form a VAR(2)
# whose true A and Sigma are the top-left N x N corners of the full ones
# (see dgp_subset()).
make_dgp_blocks <- function(spec = dgp_baseline_spec,
                            block_r2 = c(0.30, 0.60, 0.15, 0.45)) {
  A_base <- list(spec$A1, spec$A2)
  Sb     <- spec$Sigma
  nb     <- nrow(Sb)
  K      <- length(block_r2)
  N      <- nb * K
  mu_b   <- solve(diag(nb) - spec$A1 - spec$A2, spec$c_vec)

  blocks <- lapply(block_r2, function(r2) dgp_scale_to_r2(A_base, Sb, r2))

  bdiag_list <- function(mats) {
    out <- matrix(0, N, N)
    for (k in seq_along(mats)) {
      idx <- (k - 1) * nb + seq_len(nb)
      out[idx, idx] <- mats[[k]]
    }
    out
  }
  A     <- lapply(1:2, function(i) bdiag_list(lapply(blocks, `[[`, i)))
  Sigma <- bdiag_list(rep(list(Sb), K))
  mu    <- rep(mu_b, K)
  r2    <- dgp_r2(A, Sigma)

  list(
    name            = "block8",
    description     = sprintf("%d independent baseline VAR(2) blocks, block R2_trace %s",
                              K, paste(block_r2, collapse = ", ")),
    A               = A,
    p               = 2,
    Sigma           = Sigma,
    mu              = mu,
    c_vec           = as.vector((diag(N) - A[[1]] - A[[2]]) %*% mu),
    block           = rep(seq_len(K), each = nb),   # block of each series
    block_r2        = block_r2,
    block_scale     = vapply(blocks, attr, numeric(1), "scale"),
    scale           = NA_real_,
    spectral_radius = dgp_spectral_radius(A),
    R2_eq           = r2$R2_eq,
    R2_mean         = r2$R2_mean,
    R2_trace        = r2$R2_trace
  )
}

# The VAR formed by the first N series of a block DGP (N a multiple of the block size).
dgp_subset <- function(dgp, N) {
  idx <- seq_len(N)
  # whole blocks only: series N and N + 1 must lie in different blocks
  stopifnot(N <= length(dgp$block),
            N == length(dgp$block) || dgp$block[N] != dgp$block[N + 1])
  list(A = lapply(dgp$A, function(Ai) Ai[idx, idx, drop = FALSE]),
       Sigma = dgp$Sigma[idx, idx, drop = FALSE])
}

# ---- Simulation --------------------------------------------------------------

# Simulate T_obs observations after a burn-in, which removes the effect of starting at mu.
dgp_simulate <- function(dgp, T_obs, seed, burn = 500) {
  set.seed(seed)
  N     <- nrow(dgp$Sigma)
  p     <- dgp$p
  T_all <- T_obs + burn
  eps   <- matrix(rnorm(T_all * N), T_all, N) %*% chol(dgp$Sigma)
  Y     <- matrix(dgp$mu, T_all, N, byrow = TRUE)
  for (t in (p + 1):T_all) {
    y <- dgp$c_vec
    for (i in seq_len(p)) y <- y + dgp$A[[i]] %*% Y[t - i, ]
    Y[t, ] <- y + eps[t, ]
  }
  Y[(burn + 1):T_all, , drop = FALSE]
}

# ---- Standalone summary ------------------------------------------------------

if (sys.nframe() == 0L) {
  dgps <- make_dgps()
  s <- dgp_summary(dgps)
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], round, 3)
  print(s, row.names = FALSE)
  for (d in dgps) {
    cat("\n", d$name, ": ", d$description, "\n", sep = "")
    for (i in seq_len(d$p)) {
      if (any(d$A[[i]] != 0)) {
        cat(sprintf("  A%d = [%s]\n", i,
                    paste(apply(round(d$A[[i]], 3), 1, paste, collapse = " "), collapse = "; ")))
      }
    }
  }
}
