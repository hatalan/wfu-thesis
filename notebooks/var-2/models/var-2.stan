data {
    int<lower=1> T;
    int<lower=1> K;
    matrix[T,K] Y;
}

parameters {
    vector[K] c;
    matrix[K, K] A_1;
    matrix[K, K] A_2;
    vector<lower=0>[K] sigma;    // standard deviations
    corr_matrix[K] Rho;          // correlation matrix
}

model {
    c ~ normal(0, 10);
    to_vector(A_1) ~ normal(0, 10);
    to_vector(A_2) ~ normal(0, 10);

// We should talk about priors on the covariance matrix (LKJ vs. inverse-Wishart?)
    sigma ~ normal(0, 10);
    Rho ~ lkj_corr(2);
    for(t in 3:T) {
        vector[K] Y_mean = c + A_1 * Y[t-1]' + A_2 * Y[t-2]';
        Y[t]' ~ multi_normal(Y_mean, quad_form_diag(Rho, sigma));
    }
}
