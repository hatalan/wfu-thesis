functions {
  // ---- Population R2 helpers (shared by var-normal, var-minnesota, var-arr2) ----

  // Companion matrix of a VAR(p): A[i] is the lag-i coefficient matrix.
  matrix var_companion(array[] matrix A) {
    int p = size(A);
    int N = rows(A[1]);
    matrix[N * p, N * p] F = rep_matrix(0, N * p, N * p);
    for (i in 1:p)
      F[1:N, ((i - 1) * N + 1):(i * N)] = A[i];
    if (p > 1)
      F[(N + 1):(N * p), 1:(N * (p - 1))] = identity_matrix(N * (p - 1));
    return F;
  }

  // Largest eigenvalue modulus; < 1 means the VAR is stationary.
  real spectral_radius(matrix F) {
    complex_vector[rows(F)] ev = eigenvalues(F);
    real r = 0;
    for (j in 1:rows(F))
      r = fmax(r, abs(ev[j]));
    return r;
  }

  // Unconditional covariance Gamma0 = sum_j F^j Q F^j' (Q = Sigma in the top-left
  // block), solved by the doubling algorithm. Assumes spectral_radius(F) < 1.
  matrix var_gamma0(matrix F, matrix Sigma) {
    int Np = rows(F);
    int N = rows(Sigma);
    matrix[Np, Np] G = rep_matrix(0, Np, Np);
    matrix[Np, Np] Fk = F;
    G[1:N, 1:N] = Sigma;
    for (k in 1:60) {
      matrix[Np, Np] G_new = G + Fk * G * Fk';
      real change = max(abs(G_new - G));
      G = G_new;
      if (change < 1e-12 * max(abs(G)))
        break;
      Fk = Fk * Fk;
    }
    return G[1:N, 1:N];
  }

  // ---- Bayesian R2 helper (shared by var-normal, var-minnesota, var-arr2) ----

  // Sample variance, per equation, of the fitted values A_1 y_{t-1} + ... + A_p y_{t-p}
  // over the likelihood rows (p+1):T. The intercept is constant, so it drops out.
  vector var_fit_variance(matrix Y, array[] matrix A) {
    int T = rows(Y);
    int N = cols(Y);
    int p = size(A);
    matrix[T - p, N] fit = rep_matrix(0, T - p, N);
    for (i in 1:p)
      fit += Y[(p + 1 - i):(T - i)] * A[i]';
    vector[N] v;
    for (n in 1:N)
      v[n] = variance(col(fit, n));
    return v;
  }
}

data {
  int<lower=1> T;
  int<lower=1> N;                    // number of variables
  matrix[T, N] Y;
  int<lower=1> p;                    // VAR order
  real<lower=0> lambda_shape;        // Gamma hyperprior on the overall tightness lambda
  real<lower=0> lambda_rate;
  real<lower=0> theta;               // cross-variable relative tightness
  real<lower=0> lag_decay;           // prior sd shrinks as 1 / i^lag_decay
  vector[N] delta;                   // prior mean of own first lag (1 = random walk, 0 = white noise)
  vector<lower=0>[N] s;              // plug-in scale per series (e.g. univariate AR residual sd)
  real<lower=0> sigma_sd;            // sd of sigma prior
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // Prior sd of A[i][n, g] is lambda * A_scale[i][n, g]; lambda is estimated.
  array[p] matrix[N, N] A_mean;
  array[p] matrix[N, N] A_scale;
  for (i in 1:p) {
    A_mean[i] = rep_matrix(0, N, N);
    for (n in 1:N) {
      if (i == 1)
        A_mean[i][n, n] = delta[n];
      for (g in 1:N) {
        if (n == g)
          A_scale[i][n, g] = 1 / i^lag_decay;
        else
          A_scale[i][n, g] = theta * s[n] / (i^lag_decay * s[g]);
      }
    }
  }
}

parameters {
  vector[N] alpha;                   // intercept
  real<lower=0> lambda;              // overall tightness (hierarchical, as in Giannone, Lenza & Primiceri 2015)
  array[p] matrix[N, N] A_z;         // non-centered coefficients
  vector<lower=0>[N] sigma;
  corr_matrix[N] Rho;
}

transformed parameters {
  array[p] matrix[N, N] A;
  for (i in 1:p)
    A[i] = A_mean[i] + lambda * A_scale[i] .* A_z[i];
}

model {
  // priors
  target += normal_lpdf(alpha | 0, 1);
  target += gamma_lpdf(lambda | lambda_shape, lambda_rate);
  for (i in 1:p)
    target += std_normal_lpdf(to_vector(A_z[i]));
  target += normal_lpdf(sigma | 0, sigma_sd);
  target += lkj_corr_lpdf(Rho | 2);

  // likelihood
  if (!prior_only) {
    array[T - p] vector[N] y_obs;
    array[T - p] vector[N] mu;
    for (t in (p + 1):T) {
      mu[t - p] = alpha;
      for (i in 1:p)
        mu[t - p] += A[i] * Y[t - i]';
      y_obs[t - p] = Y[t]';
    }
    target += multi_normal_lpdf(y_obs | mu, quad_form_diag(Rho, sigma));
  }
}

generated quantities {
  // Population R2 implied by (A, Sigma). NaN when the draw is non-stationary.
  int<lower=0, upper=1> stable;
  vector[N] R2_eq;     // per equation: 1 - Sigma_kk / Gamma0_kk
  real R2_mean;        // equal-weight mean of R2_eq
  real R2_trace;       // 1 - tr(Sigma) / tr(Gamma0)
  {
    matrix[N, N] Sigma_gq = quad_form_diag(Rho, sigma);
    matrix[N * p, N * p] F = var_companion(A);
    stable = spectral_radius(F) < 1;
    if (stable) {
      matrix[N, N] G0 = var_gamma0(F, Sigma_gq);
      R2_eq    = 1 - diagonal(Sigma_gq) ./ diagonal(G0);
      R2_mean  = mean(R2_eq);
      R2_trace = 1 - trace(Sigma_gq) / trace(G0);
    } else {
      R2_eq    = rep_vector(not_a_number(), N);
      R2_mean  = not_a_number();
      R2_trace = not_a_number();
    }
  }

  // Bayesian R2 (Gelman, Goodrich, Gabry & Vehtari 2019, model-based residual
  // variance): fitted values on the observed lags, so it is defined for every draw,
  // stationary or not. The trace version is the quantity the ARR2 trace prior
  // targets in expectation over the coefficients.
  vector[N] R2B_eq;    // per equation: Var(fit_n) / (Var(fit_n) + sigma_n^2)
  real R2B_mean;       // equal-weight mean of R2B_eq
  real R2B_trace;      // sum_n Var(fit_n) / (sum_n Var(fit_n) + tr(Sigma))
  {
    vector[N] v  = var_fit_variance(Y, A);
    vector[N] s2 = square(sigma);
    R2B_eq    = v ./ (v + s2);
    R2B_mean  = mean(R2B_eq);
    R2B_trace = sum(v) / (sum(v) + sum(s2));
  }
}
