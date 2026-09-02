// EXNEX log-normal AFT model for right-censored data, marginalized over the
// exchangeability indicators Z_j (non-centered parameterization).
data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=1> K;
  int<lower=0> P;

  vector[N_obs] log_y_obs;
  vector[N_cens] log_y_cens;

  array[N_obs] int<lower=1, upper=K> group_obs;
  array[N_cens] int<lower=1, upper=K> group_cens;

  matrix[N_obs, P] X_obs;
  matrix[N_cens, P] X_cens;
}

parameters {
  vector[K] z;
  vector[P] beta;
  real<lower=0> sigma2;
  real<lower=0> tau2;
  real mu;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  real<lower=0> tau = sqrt(tau2);
  // Non-centered parameterization: theta = mu + tau * z with a transformed
  // EXNEX mixture prior on z so that the marginal prior on theta is unchanged.
  vector[K] theta = mu + tau * z;
}

model {
  // Priors (identical to the centered implementation and to the C++ sampler).
  sigma2 ~ inv_gamma(2.0, 2.0);
  tau2 ~ inv_gamma(2.0, 2.0);
  mu ~ normal(0.0, 100.0); // variance 1e4, SD 100
  beta ~ normal(0.0, 100.0);

  // Marginalized EXNEX mixture prior for theta_k = mu + tau * z_k.
  // The nonexchangeable component N(0, 100^2) becomes N((0 - mu)/tau, 100/tau)
  // on the z scale.
  for (k in 1:K) {
    target += log_mix(
      0.5,
      normal_lpdf(z[k] | 0.0, 1.0),
      normal_lpdf(z[k] | (0.0 - mu) / tau, 100.0 / tau)
    );
  }

  // Observed-data likelihood.
  if (P > 0) {
    log_y_obs ~ normal(theta[group_obs] + X_obs * beta, sigma);
  } else {
    log_y_obs ~ normal(theta[group_obs], sigma);
  }

  // Censored-data likelihood (Y* > Y_cens).
  if (P > 0) {
    for (i in 1:N_cens) {
      target += normal_lccdf(log_y_cens[i] | theta[group_cens[i]] + X_cens[i] * beta, sigma);
    }
  } else {
    for (i in 1:N_cens) {
      target += normal_lccdf(log_y_cens[i] | theta[group_cens[i]], sigma);
    }
  }
}
