data {
    int<lower=1> T;             // Time periods
    int<lower=1> K;             // number of lags
    real Y;
}

parameters {
    real alpha;
    vector[K] beta;
    real sigma;    // innovation
}

model {
    alpha ~ normal(0, 1);
    beta ~ normal(0, 1);
    sigma ~ normal(0, 1);

    for(t in 2:T) {
        Y[t] = alpha + beta * Y[t-1] + sigma;
    }
}
