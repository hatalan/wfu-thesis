functions {
  real arr2_ncp_lpdf(vector phi_z, vector psi, real R2, real sigma,
		 data real sigma_sd, data real mean_R2,
		 data real prec_R2, data vector cons, data real var_y) {
    return normal_lpdf(phi_z | 0, 1) +
      beta_lpdf(R2 | mean_R2 * prec_R2, (1 - mean_R2) * prec_R2) +
      normal_lpdf(sigma | 0, sigma_sd) +
      dirichlet_lpdf(psi | cons);
  }
}

data {
  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order
  // concentration vector of the Dirichlet prior
  vector<lower=0>[p] cons;
  // data for the R2D2 prior
  real<lower=0> mean_R2; // mean of the R2 prior
  real<lower=0> prec_R2; // precision of the R2 prior
  real<lower=0> sigma_sd; // sd of sigma prior
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // Variance estimate of y
  real<lower=0> var_y;
  var_y = variance(Y);
}

parameters {
  simplex[p] psi; // decomposition simplex
  real<lower=0, upper=1> R2; // coefficient of determination
  real<lower=0> sigma; // observation model sd
  real alpha; // intercept
  vector[p] phi_z;
}
transformed parameters {
  vector[p] phi;
  phi = sqrt(sigma^2 / var_y * (R2 / (1 - R2)) * psi) .* phi_z;

  vector[T] mu = rep_vector(0.0, T);
  for (t in (p+1):T) {
    mu[t] += alpha;
    for (i in 1:p) {
      mu[t] += phi[i] * Y[t-i];
    }
  }
}
model {
  // priors
  target += arr2_ncp_lpdf(phi_z | psi, R2, sigma, sigma_sd, mean_R2, prec_R2, cons, var_y);
  target += normal_lpdf(alpha | 0, 1);
  // likelihood
  if (!prior_only)
    target += normal_lpdf(Y[(p+1):T] | mu[(p+1):T], sigma);
}

generated quantities {
  vector[T]     Y_sim;
  vector[T]     mu_sim;
  
  Y_sim[1:p] = Y[1:p];
  for (t in (p + 1):T) {
    real m = alpha;
    for (i in 1:p)
      m += phi[i] * Y_sim[t - i];
    mu_sim[t] = m;
    Y_sim[t] = normal_rng(m, sigma);
  }
  
  vector[T]     Y_rep;

  Y_rep[1:p] = Y[1:p];
  for (t in (p + 1):T) {
    real m = alpha;
    for (i in 1:p)
      m += phi[i] * Y[t - i];
    Y_rep[t] = normal_rng(m, sigma);
  }
  
  real var_mu_sim = variance(mu_sim[(p+1):T]);
  real R2_sim  = var_mu_sim / (var_mu_sim + square(sigma));

  real var_mu  = variance(mu[(p+1):T]);
  real R2_data = var_mu / (var_mu + square(sigma));
}
