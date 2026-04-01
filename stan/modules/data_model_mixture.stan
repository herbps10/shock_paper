functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    real inner_sigma = epsilon_params[2];
    if(shock == 1) {
      real gamma = epsilon_params[1];
      real outer_sigma = epsilon_params[3];
      
      if(bernoulli_rng(gamma)) {
        return normal_rng(0, inner_sigma);
      }
      else {
        return normal_rng(0, outer_sigma);
      }
    }
    else {
      return normal_rng(0, inner_sigma);
    }
  }
}
transformed data {
  generate_shock_free = 1;
}
parameters {
  real<lower=0> inner_epsilon_sigma;
  real<lower=inner_epsilon_sigma> outer_epsilon_sigma;
  real<lower=0, upper=1> gamma;
}
model {
  inner_epsilon_sigma ~ normal(0, 0.75);
  outer_epsilon_sigma ~ normal(0, 10);
  gamma ~ beta(1, 50);

  target += log_mix(
    gamma,
    normal_lpdf(to_vector(diff) - to_vector(transition_function) - to_vector(shock) | 0, inner_epsilon_sigma),
    normal_lpdf(to_vector(diff) - to_vector(transition_function) - to_vector(shock) | 0, outer_epsilon_sigma)
  );
}
generated quantities {
  matrix[C, (T - 1)] likelihood_ratio;
  for(c in 1:C) {
    for(t in 1:(T - 1)) {
      likelihood_ratio[c, t] = normal_lpdf(y[c,t+1] - y[c,t] - transition_function[c,t] - shock[c,t] | 0, outer_epsilon_sigma)
        -normal_lpdf(y[c,t+1] - y[c,t] - transition_function[c,t] - shock[c,t] | 0, inner_epsilon_sigma);
    }
  }

  vector[3] epsilon_params;
  epsilon_params[1] = gamma;
  epsilon_params[2] = inner_epsilon_sigma;
  epsilon_params[3] = outer_epsilon_sigma;
}
