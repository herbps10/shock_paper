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
  //real<lower=0, upper=(shock_term == 1 ? 2 : positive_infinity())> epsilon_sigma;
  real<lower=0> epsilon_sigma;
}
model {
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);
  diff ~ normal(to_vector(transition_function) + to_vector(shock), epsilon_sigma);
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
}
