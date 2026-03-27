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
  array[hierarchical] cholesky_factor_corr[num] L_Omega_var;
}
transformed parameters {
  matrix[C, num] var;

  if(hierarchical == 1) {
    var = rep_matrix(mu_var[1]', C) + (diag_pre_multiply(sigma_var[1], L_Omega_var[1]) * raw_var')';
  }
  else {
    var = raw_var;
  }
  
  for(n in 1:num) {
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
    L_Omega_var[1] ~ lkj_corr_cholesky(1.0);
  }
  else {
    for(d in 1:D) {
      to_vector(raw_var[, d]) ~ normal(var_prior_mean[d], var_prior_sd[d]);
    }
  }
}
generated quantities {
  corr_matrix[num * hierarchical] Omega_var;
  if(hierarchical) Omega_var = multiply_lower_tri_self_transpose(L_Omega_var[1]);
}
