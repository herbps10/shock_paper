data {
  array[num] int<lower=0, upper=1> var_constrain;
  vector[num] var_lower;
  vector[num] var_upper;
  vector[num] var_prior_mean;
  vector[num] var_prior_sd;
}
parameters {
  matrix[C, num] raw_var;
  array[hierarchical] vector[num] mu_var;
  array[hierarchical] vector<lower=0>[num] sigma_var;
}
transformed parameters {
  matrix[C, num] var;

  for(n in 1:num) {
    if(hierarchical == 1) {
      var[, n] = mu_var[1][n] + sigma_var[1][n] * raw_var[, n];
    }
    else {
      var[, n] = raw_var;
    }

    if(var_constrain[n] == 1) {
      var[, n] = inv_logit(var[, n]) * (var_upper[n] - var_lower[n]) + var_lower[n];
    }
  }
}
model {
  if(hierarchical) {
    to_vector(raw_var) ~ std_normal();
    mu_var[1] ~ normal(var_prior_mean, var_prior_sd);
    sigma_var[1] ~ std_normal();
  }
  else {
    to_vector(raw_var) ~ normal(var_prior_mean, var_prior_sd);
  }
}
