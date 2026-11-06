#!/usr/bin/env Rscript
# VAR lag-misspecification sweep 02: several DGPs.
# Four coefficient priors (independent normal, Minnesota, ARR2 with flat and with
# Minnesota-decay Dirichlet concentrations) are fit at lag orders p = 1..12 to data
# from each DGP in notebooks/var/var-dgp.R (baseline VAR(2), long-decay VAR(6),
# sparse lags 1 and 4, near unit root). Fits only: prior-only runs are in var-sweep-01.
#
#   Rscript slurm/var-sweep-02.R setup       # simulate data, compile models, write grid (once)
#   Rscript slurm/var-sweep-02.R task <id>   # run one grid cell (id = SLURM_ARRAY_TASK_ID)
#   Rscript slurm/var-sweep-02.R grid        # print the grid
#   Rscript slurm/var-sweep-02.R ntasks      # print the number of array tasks
#   Rscript slurm/var-sweep-02.R dgps        # print the DGP summary
#
# Run from the repository root. Submit with slurm/submit-var-sweep-02.sh.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

cmdstan_dir <- Sys.getenv("CMDSTAN", "/deac/sta/classes/sta720/software/cmdstan/2.37.0")
if (dir.exists(cmdstan_dir)) set_cmdstan_path(cmdstan_dir)

var_dir  <- "notebooks/var"
sweep_id <- "var-sweep-02"                     # change for a new sweep version
sim_dir  <- file.path(var_dir, "simdata", sweep_id)
out_dir  <- file.path(var_dir, "results", sweep_id)

source(file.path(var_dir, "var-dgp.R"))        # make_dgps(), dgp_simulate(), dgp_summary()

# Settings
settings <- list(
  dgps        = c("baseline", "long_decay", "sparse", "near_unit_root"),
  n_rep       = 1,                            # simulated datasets per DGP
  models      = c("normal", "minnesota", "arr2_flat", "arr2_minn"),
  p_grid      = 1:12,
  modes       = c(fit = 0L),                  # value = prior_only flag; add prior = 1L for prior-only runs
  T_obs       = 200,                          # simulated length (after burn-in)
  T_train     = 150,                          # fit on rows 1:T_train, rest held out
  burn        = 500,                          # simulation burn-in, discarded
  data_seed   = 2026,                         # dataset seed = data_seed + 1000 * dgp + rep
  common_sample = TRUE,                       # same effective rows for every p
  standardize = TRUE,                         # scale with train-window mean/sd

  chains = 4, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  seed_base = 20261106,

  sigma_sd  = 1,
  normal    = list(A_sd = 1),
  minnesota = list(lambda = 0.2, theta = 0.5, lag_decay = 2,
                   delta_own = 0, s_ar_order = 4),
  # both run var-arr2.stan; they differ only in the Dirichlet concentrations
  arr2_flat = list(mean_R2 = 1/3, prec_R2 = 3, xi_split = c(1, 1),
                   conc = "flat"),
  arr2_minn = list(mean_R2 = 1/3, prec_R2 = 3, xi_split = c(1, 1),
                   conc = "minnesota")
)

### Helper functions

# Create grid from settings list (model varies fastest, DGP slowest)
make_grid <- function() {
  g <- expand.grid(model = settings$models, p = settings$p_grid,
                   mode = names(settings$modes), rep = seq_len(settings$n_rep),
                   dgp = settings$dgps, stringsAsFactors = FALSE)
  g <- g[, c("dgp", "rep", "model", "p", "mode")]
  g$task <- seq_len(nrow(g))
  g
}

sim_file <- function(dgp, rep) file.path(sim_dir, sprintf("%s_r%02d.csv", dgp, rep))

# Extract Stan model path
stan_name  <- function(model) if (startsWith(model, "arr2")) "arr2" else model
model_path <- function(model) file.path(var_dir, "models", paste0("var-", stan_name(model), ".stan"))

# Load Stan model
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

# Prepare data
prep_data <- function(dgp, rep) {
  Y <- as.matrix(read.csv(sim_file(dgp, rep), header = FALSE))
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

# Run full sweep set-up (generating data, grid)
run_setup <- function() {
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dgps <- make_dgps()[settings$dgps]
  print(dgp_summary(dgps), digits = 3, row.names = FALSE)

  # Datasets are deterministic in (data_seed, dgp, rep), so rewriting them is safe.
  for (k in seq_along(dgps)) {
    for (r in seq_len(settings$n_rep)) {
      Y <- dgp_simulate(dgps[[k]], T_obs = settings$T_obs,
                        seed = settings$data_seed + 1000 * k + r, burn = settings$burn)
      f <- sim_file(names(dgps)[k], r)
      write.table(Y, f, sep = ",", row.names = FALSE, col.names = FALSE)
      message("Wrote ", f)
    }
  }

  for (m in unique(vapply(settings$models, stan_name, ""))) load_model(m)

  grid <- make_grid()
  write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
  saveRDS(list(dgps = dgps, settings = settings), file.path(out_dir, "truth.rds"))
  message(nrow(grid), " array tasks")
}

# Run a grid cell
run_task <- function(id) {
  grid <- make_grid()
  if (id < 1 || id > nrow(grid)) stop("task id ", id, " outside 1..", nrow(grid))
  cell <- grid[id, ]
  p    <- cell$p
  prior_only <- settings$modes[[cell$mode]]

  d <- prep_data(cell$dgp, cell$rep)
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

  } else if (cell$model %in% c("arr2_flat", "arr2_minn")) {
    ar <- settings[[cell$model]]
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

  vars <- intersect(c("alpha", "A", "sigma", "Rho", "R2", "w",
                      "stable", "R2_eq", "R2_mean", "R2_trace"),
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
  f <- file.path(out_dir, sprintf("%s_r%02d_%s_p%02d_%s.rds",
                                  cell$dgp, cell$rep, cell$model, p, cell$mode))
  saveRDS(res, f)
  message("Saved ", f, " (", round(res$time, 1), "s, ",
          sum(res$diagnostics$num_divergent), " divergences)")
}

# Entry point

args <- commandArgs(trailingOnly = TRUE)
cmd  <- if (length(args)) args[1] else "usage"

switch(cmd,
  setup  = run_setup(),
  task   = run_task(as.integer(if (length(args) > 1) args[2]
                               else Sys.getenv("SLURM_ARRAY_TASK_ID"))),
  grid   = print(make_grid()),
  ntasks = cat(nrow(make_grid()), "\n"),
  dgps   = print(dgp_summary(make_dgps()[settings$dgps]), digits = 3, row.names = FALSE),
  stop("usage: Rscript slurm/var-sweep-02.R setup | task <id> | grid | ntasks | dgps")
)
