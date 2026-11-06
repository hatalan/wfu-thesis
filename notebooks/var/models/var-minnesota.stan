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
}

data {
  int<lower=1> T;
  int<lower=1> N;                    // number of variables
  matrix[T, N] Y;
  int<lower=1> p;                    // VAR order
  real<lower=0> lambda;              // overall tightness
  real<lower=0> theta;               // cross-variable relative tightness
  real<lower=0> lag_decay;           // prior sd shrinks as 1 / i^lag_decay
  vector[N] delta;                   // prior mean of own first lag (1 = random walk, 0 = white noise)
  vector<lower=0>[N] s;              // plug-in scale per series (e.g. univariate AR residual sd)
  real<lower=0> sigma_sd;            // sd of sigma prior
  int<lower=0, upper=1> prior_only;
}

transformed data {
  array[p] matrix[N, N] A_mean;
  array[p] matrix[N, N] A_sd;
  for (i in 1:p) {
    A_mean[i] = rep_matrix(0, N, N);
    for (n in 1:N) {
      if (i == 1)
        A_mean[i][n, n] = delta[n];
      for (g in 1:N) {
        if (n == g)
          A_sd[i][n, g] = lambda / i^lag_decay;
        else
          A_sd[i][n, g] = lambda * theta * s[n] / (i^lag_decay * s[g]);
      }
    }
  }
}

parameters {
  vector[N] alpha;                   // intercept
  array[p] matrix[N, N] A;
  vector<lower=0>[N] sigma;
  corr_matrix[N] Rho;
}

model {
  // priors
  target += normal_lpdf(alpha | 0, 1);
  for (i in 1:p)
    target += normal_lpdf(to_vector(A[i]) | to_vector(A_mean[i]), to_vector(A_sd[i]));
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
}
