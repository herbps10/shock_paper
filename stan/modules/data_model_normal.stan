functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
}
transformed data {
}
parameters {
  real<lower=0> epsilon_sigma;
}
model {
  epsilon_sigma ~ std_normal();
  diff ~ normal(to_vector(transition_function) + to_vector(shock), epsilon_sigma);
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
}
