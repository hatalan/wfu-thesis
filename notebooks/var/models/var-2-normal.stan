functions {
  real var_prior_lpdf(vector c, matrix A_1, matrix A_2,
                      vector sigma, matrix Rho) {
    return normal_lpdf(c | 0, 10) +
      normal_lpdf(to_vector(A_1) | 0, 10) +
      normal_lpdf(to_vector(A_2) | 0, 10) +
      normal_lpdf(sigma | 0, 10) +
      lkj_corr_lpdf(Rho | 2);
  }
}

data {
  int<lower=1> T;
  int<lower=1> K;
  matrix[T, K] Y;
}

parameters {
  vector[K] c;
  matrix[K, K] A_1;
  matrix[K, K] A_2;
  vector<lower=0>[K] sigma;    // standard deviations
  corr_matrix[K] Rho;          // correlation matrix
}

model {
  // priors
  target += var_prior_lpdf(c | A_1, A_2, sigma, Rho);

  // likelihood
  matrix[K, K] Sigma = quad_form_diag(Rho, sigma);
  for (t in 3:T) {
    vector[K] Y_mean = c + A_1 * Y[t-1]' + A_2 * Y[t-2]';
    target += multi_normal_lpdf(Y[t]' | Y_mean, Sigma);
  }
}
