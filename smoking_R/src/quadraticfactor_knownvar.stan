data {
    int<lower=0> N; // Number of samples
    int<lower=0> D; // The original dimension
    int<lower=0> K; // The latent dimension
    matrix[N, D] X_train; // The data matrix
    matrix[N, D] X_vad;
    matrix[N, D] X_all;
    matrix[N, D] holdout_mask;
    real<lower=0> data_std; // sd of data
    real<lower=0> Z_std; // sd of Z
}

parameters {
    matrix[N, K] Z; // The latent matrix
    vector[D] W0; // The weight matrix
    matrix[D, K] W1; // The weight matrix
    matrix[D, K] W2; // The weight matrix
    real<lower=0> tau; // Noise term 
    vector<lower=0>[K] alpha; // ARD prior
}

transformed parameters{
    vector<lower=0>[K] t_alpha;
    real<lower=0> t_tau;
    t_alpha = inv(sqrt(alpha));
    t_tau = inv(sqrt(tau));
}
model {
    tau ~ gamma(1,1);     
    to_vector(Z) ~ normal(0,Z_std);
    alpha ~ gamma(1e-3,1e-3);   
    W0 ~ normal(0, t_alpha[1]);
    for(k in 1:K) {
        W1[,k] ~ normal(0, t_alpha[k]);
        W2[,k] ~ normal(0, t_alpha[k]);
    }
    for(i in 1:N) to_vector(X_train[i]) ~ normal(to_vector((W0' + Z[i]*W1' + (Z[i].*Z[i])*W2').*(1-holdout_mask[i])), data_std);
} 



generated quantities{
    vector[N] rep_lp;
    vector[N] vad_lp;
    vector[N] ps;
    real predmean;
    matrix[N, D] X_pred;

    for(i in 1:N) {
        ps[i] = normal_lpdf(to_vector(X_all[i]) | to_vector(W0' + Z[i]*W1' + (Z[i].*Z[i])*W2'), data_std);
        rep_lp[i] = 0;
        vad_lp[i] = 0;   
        for(j in 1:D){
            predmean = W0[j] + Z[i,]*to_vector(W1[j,]) + (Z[i,].*Z[i,])*to_vector(W2[j,]);
            X_pred[i,j] = normal_rng(predmean, data_std);
            rep_lp[i] += normal_lpdf(X_pred[i,j] | predmean, data_std) * holdout_mask[i,j];
            vad_lp[i] += normal_lpdf(X_vad[i,j] | predmean, data_std) * holdout_mask[i,j];
        }
    }
}


