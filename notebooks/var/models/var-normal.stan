data {
  int<lower=1> T;
  int<lower=1> N;                    // number of variables
  matrix[T, N] Y;
  int<lower=1> p;                    // VAR order
  real<lower=0> A_sd;                // sd of the independent normal prior on A
  real<lower=0> sigma_sd;            // sd of sigma prior
  int<lower=0, upper=1> prior_only;
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
    target += normal_lpdf(to_vector(A[i]) | 0, A_sd);
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
