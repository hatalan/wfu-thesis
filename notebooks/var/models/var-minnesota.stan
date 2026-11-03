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
