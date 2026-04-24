data {
    int<lower=1> T;
    int<lower=1> K;
    int<lower=1> P;
    matrix[T, K] Y;
}

transformed data {
    vector[K] delta = rep_vector(1.0, K);
    real alpha_lag = 2.0;
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
    lambda ~ gamma(2, 10);
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
                real prior_var;
                
                if (i == j && ell == 1) {
                    prior_mean = 1; // defining diagonal mean priors
                } else {
                    prior_mean = 0; // defining off-diagonal mean priors
                  }
                
                if (i == j) {
                    prior_var = lambda^2 / ell^alpha_lag * (1 / sigma[i]^2); // defining diagonal variance priors
                } else {
                    prior_var = lambda^2 * theta^2 / ell^alpha_lag * (sigma[j]^2 / sigma[i]^2); // defining off-diagonal variance priors
                }
                
                A[ell][i, j] ~ normal(prior_mean, sqrt(prior_var));
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