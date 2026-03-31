functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
}
data {
  real<lower=0> epsilon_sigma_prior_mu;
  real<lower=0> epsilon_sigma_prior_sd;
}
transformed data {
}
parameters {
  real<lower=0> epsilon_sigma;
}
model {
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);
  to_vector(y[, 2:T] - y[, 1:(T - 1)]) ~ normal(to_vector(transition_function + shock[, 2:T] - shock[, 1:(T - 1)]), epsilon_sigma);
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
}
