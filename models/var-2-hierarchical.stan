// Hierarchical VAR(P) with a Minnesota prior.
//
// Cluster/sweep version of notebooks/var-2-hierarchical.stan: the gamma
// hyperprior on the overall tightness `lambda` is parameterized via data
// (lambda_shape, lambda_rate) so a SLURM array can sweep prior specifications.
// Setting lambda_shape=2, lambda_rate=10 reproduces the notebook's gamma(2, 10).

data {
    int<lower=1> T;
    int<lower=1> K;
    int<lower=1> P;
    matrix[T, K] Y;

    // Gamma hyperprior on lambda (the swept quantity): lambda ~ gamma(shape, rate)
    real<lower=0> lambda_shape;
    real<lower=0> lambda_rate;
}

transformed data {
    vector[K] delta = rep_vector(1.0, K);
    real alpha_lag = 1.0;
}

parameters {
    vector[K] c;

    array[P] matrix[K, K] A;

    vector<lower=0>[K] sigma;
    corr_matrix[K] Rho;

    real<lower=0> lambda;        // overall tightness
    real<lower=0, upper=1> theta; // cross-equation relative tightness
}

model {
    // --- Hyperpriors ---
    lambda ~ gamma(lambda_shape, lambda_rate);
    theta ~  beta(2, 2);

    // --- Priors on innovation scale and correlation ---
    sigma ~ normal(0, 1);
    Rho ~ lkj_corr(2);

    // Prior on intercept
    c ~ normal(0, 10);

    // Minnesota prior
    for (ell in 1:P) {
        for (i in 1:K) {
            for (j in 1:K) {
                real prior_mean;
                real prior_sd;

                if (i == j && ell == 1) {
                    prior_mean = 1; // defining diagonal mean priors
                } else {
                    prior_mean = 0; // defining off-diagonal mean priors
                  }

                if (i == j) {
                    prior_sd = lambda / ell^alpha_lag; // defining diagonal variance priors
                } else {
                    prior_sd = (lambda * theta * sigma[i]) / (ell^alpha_lag * sigma[j]); // defining off-diagonal variance priors
                }

                A[ell][i, j] ~ normal(prior_mean, prior_sd);
            }
        }
    }

    for (t in (P+1):T) {
        vector[K] Y_mean = c;

        for (ell in 1:P) {
            Y_mean += A[ell] * Y[t - ell]';
        }

        Y[t]' ~ multi_normal(Y_mean, quad_form_diag(Rho, sigma));
    }
}
