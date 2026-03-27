data {
  int<lower=0, upper=1> var_constrain;
  real var_lower;
  real var_upper;
  real var_sigma_lower;
  real var_prior_mean;
  real var_prior_sd;
}
parameters {
  vector[C] raw_var;
  array[hierarchical] real mu_var;
  array[hierarchical] real sigma_var;
}
transformed parameters {
  vector[C] var;
  
  if(hierarchical == 1) {
    var = mu_var[1] + (var_sigma_lower + exp(sigma_var[1])) * raw_var;
  }
  else {
    var = raw_var;
  }
  if(var_constrain == 1) {
    var = inv_logit(var) * (var_upper - var_lower) + var_lower;
  }
}
model {
  if(hierarchical) {
    raw_var ~ std_normal();
    mu_var[1] ~ normal(var_prior_mean, var_prior_sd);
    //sigma_var[1] ~ std_normal();
    sigma_var[1] ~ normal(-1, 0.5);
  }
  else {
    raw_var ~ normal(var_prior_mean, var_prior_sd);
  }
}
