data {
  int<lower=0> N;
  array[N] int<lower=1, upper=3> income;
  array[N] int<lower=1, upper=3> education;
  array[N] int<lower=1, upper=3> age;
  array[N] int<lower=1, upper=51> state;  
  array[N] int<lower=0> y;
  array[3, 3, 3, 51] real<lower=0> P;
}

parameters {
  real alpha;
  real<lower=0> sigma_beta_1;
  vector<multiplier=sigma_beta_1>[3] beta_1;
  real<lower=0> sigma_beta_2;
  vector<multiplier=sigma_beta_2>[3] beta_2;
  real<lower=0> sigma_beta_3;
  vector<multiplier=sigma_beta_3>[3] beta_3;
  real<lower=0> sigma_beta_4;
  vector<multiplier=sigma_beta_4>[51] beta_4;
}

model {
  y ~ bernoulli_logit(alpha + beta_1[income] + beta_2[education] + beta_3[age] + beta_4[state]);
  alpha ~ normal(0, 2);
  beta_1 ~ normal(0, sigma_beta_1);
  beta_2 ~ normal(0, sigma_beta_2);
  beta_3 ~ normal(0, sigma_beta_3);
  beta_4 ~ normal(0, sigma_beta_4);
  { sigma_beta_1, sigma_beta_2, sigma_beta_3, sigma_beta_4 } ~ normal(0, 1);
}

generated quantities {
  real expect_pos = 0;
  real total = 0;
  for (b in 1:3) {
    for (c in 1:3) {
      for (d in 1:3) {
        for( e in 1:51) {
          total += P[b, c, d, e];
          expect_pos += P[b, c, d]
            * inv_logit(alpha + beta_1[b] + beta_2[c] + beta_3[d] + beta_4[e]);
        }
      }
    }
  }
  real<lower=0, upper=1> phi = expect_pos / total;
}
