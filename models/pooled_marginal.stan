// Complete-pooling log-normal AFT model for right-censored data.
data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=1> K;
  int<lower=0> P;

  vector[N_obs] log_y_obs;
  vector[N_cens] log_y_cens;

  matrix[N_obs, P] X_obs;
  matrix[N_cens, P] X_cens;
}

parameters {
  real theta; // Complete pooling: a single shared intercept.
  vector[P] beta;
  real<lower=0> sigma2;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
}

model {
  sigma2 ~ inv_gamma(2.0, 2.0);
  theta ~ normal(0.0, 100.0);
  beta ~ normal(0.0, 100.0);

  // Observed-data likelihood.
  if (P > 0) {
    log_y_obs ~ normal(theta + X_obs * beta, sigma);
  } else {
    log_y_obs ~ normal(theta, sigma);
  }

  // Censored-data likelihood.
  if (P > 0) {
    for (i in 1:N_cens) {
      target += normal_lccdf(log_y_cens[i] | theta + X_cens[i] * beta, sigma);
    }
  } else {
    for (i in 1:N_cens) {
      target += normal_lccdf(log_y_cens[i] | theta, sigma);
    }
  }
}
