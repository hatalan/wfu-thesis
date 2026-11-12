#!/usr/bin/env Rscript
# VAR sweep 03: number of variables.
# Four coefficient priors (independent normal, hierarchical Minnesota, ARR2 with flat and
# with Minnesota-decay Dirichlet concentrations) are fit to the first N = 2, 4, 6, 8 series
# of an 8-variable VAR(2) made of four independent 2-variable blocks (make_dgp_blocks() in
# notebooks/var/var-dgp.R), at p = 1, 2, 4, with and without the likelihood, on 3 simulated
# datasets. Each subset is made of whole blocks, so it is itself a VAR(2) with a known truth.
#
#   Rscript slurm/var-sweep-03.R setup       # simulate data, compile models, write grid (once)
#   Rscript slurm/var-sweep-03.R task <id>   # run one grid cell (id = SLURM_ARRAY_TASK_ID)
#   Rscript slurm/var-sweep-03.R grid        # print the grid
#   Rscript slurm/var-sweep-03.R ntasks      # print the number of array tasks
#   Rscript slurm/var-sweep-03.R dgps        # print the DGP summary
#
# Run from the repository root. Submit with slurm/submit-var-sweep-03.sh.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

cmdstan_dir <- Sys.getenv("CMDSTAN", "/deac/sta/classes/sta720/software/cmdstan/2.37.0")
if (dir.exists(cmdstan_dir)) set_cmdstan_path(cmdstan_dir)

var_dir  <- "notebooks/var"
sweep_id <- "var-sweep-03"                     # change for a new sweep version
sim_dir  <- file.path(var_dir, "simdata", sweep_id)
out_dir  <- file.path(var_dir, "results", sweep_id)

source(file.path(var_dir, "var-dgp.R"))        # make_dgp_blocks(), dgp_subset(), dgp_simulate()

# Settings
settings <- list(
  block_r2    = c(0.30, 0.60, 0.15, 0.45),    # population R2_trace of blocks 1-4
  N_grid      = c(2, 4, 6, 8),                # fit the first N series (whole blocks)
  n_rep       = 3,                            # simulated datasets
  models      = c("normal", "minnesota", "arr2_flat", "arr2_minn"),
  p_grid      = c(1, 2, 4),
  modes       = c(fit = 0L, prior = 1L),      # value = prior_only flag
  T_obs       = 200,                          # simulated length (after burn-in)
  T_train     = 150,                          # fit on rows 1:T_train
  burn        = 500,                          # simulation burn-in, discarded
  data_seed   = 3026,                         # dataset seed = data_seed + rep
  common_sample = TRUE,                       # same effective rows for every p
  standardize = TRUE,                         # scale with train-window mean/sd

  chains = 4, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  seed_base = 20261112,

  sigma_sd  = 1,
  normal    = list(A_sd = 1),
  # Hierarchical Minnesota: lambda ~ Gamma with mode 0.2 and sd 0.4 (Giannone, Lenza &
  # Primiceri 2015); theta ~ Beta(2, 2), centred on the conventional 0.5; lag_decay ~ Gamma
  # with mode 1 (harmonic decay, the usual default) and sd 0.5, covering the usual 1-2 range.
  minnesota = list(lambda_mode = 0.2, lambda_sd = 0.4, theta_a = 2, theta_b = 2,
                   lag_decay_mode = 1, lag_decay_sd = 0.5,
                   delta_own = 0, s_ar_order = 4),
  # both run var-arr2.stan; they differ only in the Dirichlet concentrations
  arr2_flat = list(mean_R2 = 1/3, prec_R2 = 3, xi_split = c(1, 1),
                   conc = "flat"),
  arr2_minn = list(mean_R2 = 1/3, prec_R2 = 3, xi_split = c(1, 1),
                   conc = "minnesota")
)

### Helper functions

# Gamma shape and rate from a mode and standard deviation (mode > 0)
gamma_mode_sd <- function(mode, sd) {
  s <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd)   # sqrt(shape)
  c(shape = s^2, rate = s / sd)
}

make_dgp <- function() make_dgp_blocks(block_r2 = settings$block_r2)

# Create grid from settings list (model varies fastest, replicate slowest)
make_grid <- function() {
  g <- expand.grid(model = settings$models, p = settings$p_grid, N = settings$N_grid,
                   mode = names(settings$modes), rep = seq_len(settings$n_rep),
                   stringsAsFactors = FALSE)
  g <- g[, c("rep", "mode", "N", "p", "model")]
  g$task <- seq_len(nrow(g))
  g
}

sim_file <- function(rep) file.path(sim_dir, sprintf("block8_r%02d.csv", rep))

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

# Prepare data: the first N series, standardized with the train window
prep_data <- function(rep, N) {
  Y <- as.matrix(read.csv(sim_file(rep), header = FALSE))[, seq_len(N), drop = FALSE]
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

dgp_table <- function(dgp) {
  data.frame(series = seq_along(dgp$block), block = dgp$block,
             block_R2_trace = dgp$block_r2[dgp$block], R2_eq = dgp$R2_eq)
}

# Run full sweep set-up (generating data, grid)
run_setup <- function() {
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dgp <- make_dgp()
  print(dgp_table(dgp), digits = 3, row.names = FALSE)
  message("spectral radius ", round(dgp$spectral_radius, 3))

  # Datasets are deterministic in (data_seed, rep), so rewriting them is safe.
  for (r in seq_len(settings$n_rep)) {
    Y <- dgp_simulate(dgp, T_obs = settings$T_obs,
                      seed = settings$data_seed + r, burn = settings$burn)
    write.table(Y, sim_file(r), sep = ",", row.names = FALSE, col.names = FALSE)
    message("Wrote ", sim_file(r))
  }

  for (m in unique(vapply(settings$models, stan_name, ""))) load_model(m)

  grid <- make_grid()
  write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)
  saveRDS(list(dgp = dgp, settings = settings), file.path(out_dir, "truth.rds"))
  message(nrow(grid), " array tasks")
}

# Run a grid cell
run_task <- function(id) {
  grid <- make_grid()
  if (id < 1 || id > nrow(grid)) stop("task id ", id, " outside 1..", nrow(grid))
  cell <- grid[id, ]
  p    <- cell$p
  N    <- cell$N
  prior_only <- settings$modes[[cell$mode]]

  d <- prep_data(cell$rep, N)
  T_train <- settings$T_train
  p_max   <- max(settings$p_grid)
  first   <- if (settings$common_sample) p_max - p + 1 else 1
  Y_fit   <- d$Y[first:T_train, , drop = FALSE]

  stan_data <- list(T = nrow(Y_fit), N = N, Y = Y_fit, p = p,
                    sigma_sd = settings$sigma_sd, prior_only = prior_only)

  if (cell$model == "normal") {
    stan_data$A_sd <- settings$normal$A_sd

  } else if (cell$model == "minnesota") {
    mn <- settings$minnesota
    Y_tr <- d$Y[seq_len(T_train), , drop = FALSE]
    s <- apply(Y_tr, 2, function(y)
      sqrt(arima(y, order = c(mn$s_ar_order, 0, 0))$sigma2))
    lam <- gamma_mode_sd(mn$lambda_mode, mn$lambda_sd)
    dec <- gamma_mode_sd(mn$lag_decay_mode, mn$lag_decay_sd)
    stan_data <- c(stan_data, list(lambda_shape    = lam[["shape"]],
                                   lambda_rate     = lam[["rate"]],
                                   theta_a         = mn$theta_a,
                                   theta_b         = mn$theta_b,
                                   lag_decay_shape = dec[["shape"]],
                                   lag_decay_rate  = dec[["rate"]],
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

  vars <- intersect(c("alpha", "A", "sigma", "Rho", "R2", "w", "lambda", "theta", "lag_decay",
                      "stable", "R2_eq", "R2_mean", "R2_trace",
                      "R2B_eq", "R2B_mean", "R2B_trace"),
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
  f <- file.path(out_dir, sprintf("block8_r%02d_N%02d_%s_p%02d_%s.rds",
                                  cell$rep, N, cell$model, p, cell$mode))
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
  dgps   = {
    dgp <- make_dgp()
    print(dgp_table(dgp), digits = 3, row.names = FALSE)
    cat("spectral radius", round(dgp$spectral_radius, 3), "\n")
    for (N in settings$N_grid) {
      s <- dgp_subset(dgp, N)
      cat(sprintf("N = %d: population R2_trace %.3f\n", N, dgp_r2(s$A, s$Sigma)$R2_trace))
    }
  },
  stop("usage: Rscript slurm/var-sweep-03.R setup | task <id> | grid | ntasks | dgps")
)
