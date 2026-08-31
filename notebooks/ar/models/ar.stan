data {
  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=1> p; // AR order
  real<lower=0> sigma_sd; // sd of sigma prior (shared with arr2.stan)
  real<lower=0> phi_sd; // sd of the naive independent prior on phi
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
  target += normal_lpdf(Y[(p+1):T] | mu[(p+1):T], sigma);
}

generated quantities {
  vector[T - p] log_lik;
  vector[T - p] Y_rep;
  vector[T]     Y_sim;

  for (t in (p + 1):T) {
    log_lik[t - p] = normal_lpdf(Y[t] | mu[t], sigma);
    Y_rep[t - p]   = normal_rng(mu[t], sigma);
  }

  Y_sim[1:p] = Y[1:p];
  for (t in (p + 1):T) {
    real m = alpha;
    for (i in 1:p)
      m += phi[i] * Y_sim[t - i];
    Y_sim[t] = normal_rng(m, sigma);
  }
}
