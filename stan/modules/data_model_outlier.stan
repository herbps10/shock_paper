functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
}
data {
  real<lower=0> outlier_threshold;
}
transformed data {
  int n_below_threshold = 0;
  for(i in 1:(C * (T - 1))) {
    n_below_threshold += (abs(diff[i]) < outlier_threshold) ? 1 : 0;
  }
  array[n_below_threshold] int indices_below_threshold;
  
  {
    int index = 1;
    for(i in 1:(C * (T - 1))) {
      if(abs(diff[i]) < outlier_threshold) {
        indices_below_threshold[index] = i;
        index += 1;
      }
    }
  }
}
parameters {
  real<lower=0> epsilon_sigma;
}
model {
  epsilon_sigma ~ std_normal();

  if(outlier_threshold < 1000) {
    diff[indices_below_threshold] ~ normal(to_vector(transition_function)[indices_below_threshold], epsilon_sigma);
  }
  else {
    diff ~ normal(to_vector(transition_function), epsilon_sigma);
  }
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
}
