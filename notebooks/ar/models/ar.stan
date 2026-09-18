data {
  int<lower=1> T; // length of time series
  vector[T] Y; // observations
  int<lower=1> p; // AR order
  real<lower=0> sigma_sd; // sd of sigma prior (shared with arr2.stan)
  real<lower=0> phi_sd; // sd of the naive independent prior on phi
  int<lower=0, upper=1> prior_only;
}

parameters {
  real alpha; // intercept
  vector[p] phi; // AR coefficients
  real<lower=0> sigma; // observation model sd
}

transformed parameters {
  vector[T] mu = rep_vector(0.0, T);
  for (t in (p+1):T) {
    mu[t] = alpha;
    for (i in 1:p) {
      mu[t] += phi[i] * Y[t-i];
    }
  }
}

model {
  // priors
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(phi | 0, phi_sd);
  target += normal_lpdf(sigma | 0, sigma_sd);
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
