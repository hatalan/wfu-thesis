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

  real arr2_var_lpdf(vector phi_z, vector beta_z,
                     vector w, vector psi_raw, vector gamma_raw,
                     real R2, vector sigma, matrix Rho,
                     data real sigma_sd, data real mean_R2, data real prec_R2,
                     data vector xi_split, data vector xi_own, data vector xi_cross) {
    return normal_lpdf(phi_z | 0, 1) +
      normal_lpdf(beta_z | 0, 1) +
      beta_lpdf(R2 | mean_R2 * prec_R2, (1 - mean_R2) * prec_R2) +
      normal_lpdf(sigma | 0, sigma_sd) +
      lkj_corr_lpdf(Rho | 2) +
      dirichlet_lpdf(w | xi_split) +
      dirichlet_lpdf(psi_raw | xi_own) +
      dirichlet_lpdf(gamma_raw | xi_cross);
  }
}

data {
  int<lower=1> T;
  int<lower=1> N; // number of variables
  matrix[T, N] Y;
  int<lower=1> p; // VAR order
  // concentration vectors of the Dirichlet priors
  vector<lower=0>[2] xi_split;
  vector<lower=0>[N * p] xi_own;
  vector<lower=0>[N * (N - 1) * p] xi_cross;
  // data for the R2 prior
  real<lower=0> mean_R2; // mean of the R2 prior
  real<lower=0> prec_R2; // precision of the R2 prior
  real<lower=0> sigma_sd; // sd of sigma prior
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // plug-in variance of each series (diagonal of the sample covariance S)
  vector[N] diag_S;
  for (n in 1:N) diag_S[n] = variance(col(Y, n));

  // plug-in variance of each coefficient's regressor
  vector[N * p] S_phi;
  vector[N * (N - 1) * p] S_beta;
  for (n in 1:N) {
    for (i in 1:p) {
      S_phi[(n - 1) * p + i] = diag_S[n];
      for (g in 1:N) {
        if (g != n) {
          int gg = g < n ? g : g - 1;
          S_beta[(n - 1) * (N - 1) * p + (i - 1) * (N - 1) + gg] = diag_S[g];
        }
      }
    }
  }
}

parameters {
  simplex[2] w;
  simplex[N * p] psi_raw;
  simplex[N * (N - 1) * p] gamma_raw;
  real<lower=0, upper=1> R2; // coefficient of determination
  vector<lower=0>[N] sigma;
  corr_matrix[N] Rho;
  vector[N] alpha; // intercept
  vector[N * p] phi_z;
  vector[N * (N - 1) * p] beta_z;
}

transformed parameters {
  real tau2     = R2 / (1 - R2);
  real tr_Sigma = dot_self(sigma);

  vector[N * p] psi = w[1] * psi_raw;
  vector[N * (N - 1) * p] gamma = w[2] * gamma_raw;

  vector[N * p] phi = sqrt(tr_Sigma * tau2 * psi   ./ S_phi)  .* phi_z;
  vector[N * (N - 1) * p] beta = sqrt(tr_Sigma * tau2 * gamma ./ S_beta) .* beta_z;

  array[p] matrix[N, N] A;
  for (i in 1:p) {
    for (n in 1:N) {
      A[i][n, n] = phi[(n - 1) * p + i];
      for (g in 1:N) {
        if (g != n) {
          int gg = g < n ? g : g - 1;
          A[i][n, g] = beta[(n - 1) * (N - 1) * p + (i - 1) * (N - 1) + gg];
        }
      }
    }
  }

  matrix[T, N] mu = rep_matrix(0.0, T, N);
  for (t in (p + 1):T) {
    mu[t] = alpha';
    for (i in 1:p) {
      mu[t] += Y[t - i] * A[i]';
    }
  }
}

model {
  // priors
  target += arr2_var_lpdf(phi_z | beta_z, w, psi_raw, gamma_raw, R2, sigma, Rho,
                          sigma_sd, mean_R2, prec_R2, xi_split, xi_own, xi_cross);
  target += normal_lpdf(alpha | 0, 1);

  // likelihood
  if (!prior_only) {
    matrix[N, N] Sigma = quad_form_diag(Rho, sigma);
    for (t in (p + 1):T)
      target += multi_normal_lpdf(Y[t] | mu[t], Sigma);
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
