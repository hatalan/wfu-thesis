#!/usr/bin/env Rscript
# VAR lag-misspecification sweep.
# Three coefficient priors (independent normal, Minnesota, ARR2) are fit to one
# simulated VAR(2) dataset at lag orders p = 1..8, both to the data and prior-only.
#
#   Rscript slurm/var-sweep.R setup       # simulate data, compile models, write grid (once)
#   Rscript slurm/var-sweep.R task <id>   # run one grid cell (id = SLURM_ARRAY_TASK_ID)
#   Rscript slurm/var-sweep.R grid        # print the grid
#
# Run from the repository root. Submit with slurm/submit-var-sweep.sh.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

cmdstan_dir <- Sys.getenv("CMDSTAN", "/deac/sta/classes/sta720/software/cmdstan/2.37.0")
if (dir.exists(cmdstan_dir)) set_cmdstan_path(cmdstan_dir)

var_dir  <- "notebooks/var"
sim_file <- file.path(var_dir, "simdata", "Y_sim.csv")
out_dir  <- file.path(var_dir, "results", "var-sweep")

# ---- Settings (decisions live here) -----------------------------------------

settings <- list(
  models      = c("normal", "minnesota", "arr2"),
  p_grid      = 1:8,
  modes       = c(fit = 0L, prior = 1L),      # value = prior_only flag
  T_train     = 150,                          # fit on rows 1:T_train, rest held out
  common_sample = TRUE,                       # same effective rows for every p
  standardize = TRUE,                         # scale with train-window mean/sd

  chains = 4, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  seed_base = 20261103,

  sigma_sd  = 1,                              # shared by all three models
  normal    = list(A_sd = 1),
  minnesota = list(lambda = 0.2, theta = 0.5, lag_decay = 2,
                   delta_own = 0, s_ar_order = 4),
  arr2      = list(mean_R2 = 1/3, prec_R2 = 3, xi_split = c(1, 1),
                   conc = "flat")             # "flat" or "minnesota"
)

# DGP: identical to var-2-notebook.qmd, so Y_sim.csv is the same file either way.
dgp <- list(
  T_obs = 200, seed = 1,
  c_vec = c(2.5, 1.5),
  A = list(matrix(c(0.45, -0.15, -0.35, 0.30), nrow = 2),
           matrix(c(0.05, -0.05, -0.05, 0.03), nrow = 2)),
  Sigma = matrix(c(1.0, 0.1, 0.1, 1.0), nrow = 2)
)

# ---- Helpers ----------------------------------------------------------------

make_grid <- function() {
  g <- expand.grid(model = settings$models, p = settings$p_grid,
                   mode = names(settings$modes), stringsAsFactors = FALSE)
  g$task <- seq_len(nrow(g))
  g
}

sim_var <- function(c_vec, A_1, A_2, Sigma, T_obs, seed) {
  set.seed(seed)
  K <- length(c_vec)
  Y <- matrix(0, nrow = T_obs, ncol = K)
  epsilon <- MASS::mvrnorm(n = T_obs, mu = rep(0, K), Sigma = Sigma)
  mu <- solve(diag(K) - A_1 - A_2, c_vec)
  Y[1, ] <- MASS::mvrnorm(1, mu, Sigma)
  Y[2, ] <- MASS::mvrnorm(1, mu, Sigma)
  for (t in 3:T_obs)
    Y[t, ] <- c_vec + A_1 %*% Y[t - 1, ] + A_2 %*% Y[t - 2, ] + epsilon[t, ]
  Y
}

model_path <- function(model) file.path(var_dir, "models", paste0("var-", model, ".stan"))

load_model <- function(model) {
  cmdstan_model(model_path(model), stanc_options = list("O1"))  # no-op if compiled
}

# ARR2 Dirichlet concentrations, laid out to match the indexing in var-arr2.stan
arr2_conc <- function(p, N, type) {
  base <- switch(type,
                 flat      = rep(1, p),
                 minnesota = (p^2 / 10) / (1:p)^2)
  list(xi_own   = rep(base, N),
       xi_cross = rep(rep(base, each = N - 1), N))
}

prep_data <- function() {
  Y <- as.matrix(read.csv(sim_file, header = FALSE))
  Y_tr <- Y[seq_len(settings$T_train), , drop = FALSE]
  if (settings$standardize) {
    center <- colMeans(Y_tr)
    scale  <- apply(Y_tr, 2, sd)
  } else {
    center <- rep(0, ncol(Y)); scale <- rep(1, ncol(Y))
  }
  Ys <- sweep(sweep(Y, 2, center), 2, scale, "/")
  list(Y = Ys, center = center, scale = scale)
}

# ---- Setup: simulate + compile ----------------------------------------------

run_setup <- function() {
  dir.create(dirname(sim_file), showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(sim_file)) {
    Y <- sim_var(dgp$c_vec, dgp$A[[1]], dgp$A[[2]], dgp$Sigma, dgp$T_obs, dgp$seed)
    write.table(Y, sim_file, sep = ",", row.names = FALSE, col.names = FALSE)
    message("Wrote ", sim_file)
  } else {
    message("Using existing ", sim_file)
  }

  for (m in settings$models) load_model(m)

  grid <- make_grid()
  write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
  saveRDS(list(dgp = dgp, settings = settings), file.path(out_dir, "truth.rds"))
  message(nrow(grid), " array tasks")
}

# ---- One grid cell ----------------------------------------------------------

run_task <- function(id) {
  grid <- make_grid()
  if (id < 1 || id > nrow(grid)) stop("task id ", id, " outside 1..", nrow(grid))
  cell <- grid[id, ]
  p    <- cell$p
  prior_only <- settings$modes[[cell$mode]]

  d <- prep_data()
  T_train <- settings$T_train
  p_max   <- max(settings$p_grid)
  first   <- if (settings$common_sample) p_max - p + 1 else 1
  Y_fit   <- d$Y[first:T_train, , drop = FALSE]
  N       <- ncol(Y_fit)

  stan_data <- list(T = nrow(Y_fit), N = N, Y = Y_fit, p = p,
                    sigma_sd = settings$sigma_sd, prior_only = prior_only)

  s <- NULL
  if (cell$model == "normal") {
    stan_data$A_sd <- settings$normal$A_sd

  } else if (cell$model == "minnesota") {
    mn <- settings$minnesota
    Y_tr <- d$Y[seq_len(T_train), , drop = FALSE]
    s <- apply(Y_tr, 2, function(y)
      sqrt(arima(y, order = c(mn$s_ar_order, 0, 0))$sigma2))
    stan_data <- c(stan_data, list(lambda = mn$lambda, theta = mn$theta,
                                   lag_decay = mn$lag_decay,
                                   delta = rep(mn$delta_own, N), s = s))

  } else if (cell$model == "arr2") {
    ar <- settings$arr2
    stan_data <- c(stan_data,
                   list(mean_R2 = ar$mean_R2, prec_R2 = ar$prec_R2,
                        xi_split = ar$xi_split),
                   arr2_conc(p, N, ar$conc))
  }

  mod <- load_model(cell$model)
  fit <- mod$sample(
    data            = stan_data,
    chains          = settings$chains,
    parallel_chains = settings$chains,
    iter_warmup     = settings$iter_warmup,
    iter_sampling   = settings$iter_sampling,
    adapt_delta     = settings$adapt_delta,
    seed            = settings$seed_base + id,
    refresh         = 0,
    show_messages   = FALSE,
    show_exceptions = FALSE
  )

  vars <- intersect(c("alpha", "A", "sigma", "Rho", "R2", "w"),
                    fit$metadata()$stan_variables)

  res <- list(
    cell        = cell,
    stan_data   = stan_data[setdiff(names(stan_data), "Y")],
    rows_fit    = c(first, T_train),
    center      = d$center,
    scale       = d$scale,
    draws       = fit$draws(vars, format = "draws_df"),
    summary     = fit$summary(vars),
    diagnostics = fit$diagnostic_summary(quiet = TRUE),
    time        = fit$time()$total
  )

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(out_dir, sprintf("%s_p%02d_%s.rds", cell$model, p, cell$mode))
  saveRDS(res, f)
  message("Saved ", f, " (", round(res$time, 1), "s, ",
          sum(res$diagnostics$num_divergent), " divergences)")
}

# ---- Entry point ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
cmd  <- if (length(args)) args[1] else "usage"

switch(cmd,
  setup = run_setup(),
  task  = run_task(as.integer(if (length(args) > 1) args[2]
                              else Sys.getenv("SLURM_ARRAY_TASK_ID"))),
  grid  = print(make_grid()),
  stop("usage: Rscript slurm/var-sweep.R setup | task <id> | grid")
)
