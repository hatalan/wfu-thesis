functions {
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
  matrix[N, N] Sigma = quad_form_diag(Rho, sigma);
  for (t in (p + 1):T)
    target += multi_normal_lpdf(Y[t] | mu[t], Sigma);
}
