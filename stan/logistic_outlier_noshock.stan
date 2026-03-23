functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
  
  vector rate_double_logistic(vector x, vector Delta1, vector Delta2,
                              vector Delta3, vector Delta4, vector k,
                              vector z) {
    real A1 = 4.4;
    real A2 = 0.5;
    return k
           .* inv(1 + exp(-A1 .* inv(Delta2) .* (x - Delta1 - A2 * Delta2)))
           + (z - k)
             .* inv(
                    1
                    + exp(
                          -A1 * inv(Delta4)
                          .* (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
  }
  
  vector rep_vector_times(vector x, int times) {
    int s = size(x);
    vector[s * times] y;
    for (i in 1 : times) {
      y[((i - 1) * s + 1) : (i * s)] = x;
    }
    return y;
  }
}
data {
  int C;
  int T;
  int Tpred;
  
  matrix[C, T] y;
  
  int num_grid;
  vector[num_grid] grid;
  
  int<lower=0, upper=1> include_prior;
  int<lower=0, upper=1> hierarchical;
  
  real<lower=0> outlier_threshold;
  
  int<lower=0, upper=1> Delta1_constrain;
  real Delta1_lower;
  real Delta1_upper;
  real Delta1_prior_mean;
  real Delta1_prior_sd;
  int<lower=0, upper=1> Delta2_constrain;
  real Delta2_lower;
  real Delta2_upper;
  real Delta2_prior_mean;
  real Delta2_prior_sd;
  int<lower=0, upper=1> Delta3_constrain;
  real Delta3_lower;
  real Delta3_upper;
  real Delta3_prior_mean;
  real Delta3_prior_sd;
  int<lower=0, upper=1> Delta4_constrain;
  real Delta4_lower;
  real Delta4_upper;
  real Delta4_prior_mean;
  real Delta4_prior_sd;
  int<lower=0, upper=1> k_constrain;
  real k_lower;
  real k_upper;
  real k_prior_mean;
  real k_prior_sd;
  int<lower=0, upper=1> z_constrain;
  real z_lower;
  real z_upper;
  real z_prior_mean;
  real z_prior_sd;
}
transformed data {
  vector[C * (T - 1)] diff = to_vector(y[ : , 2 : T] - y[ : , 1 : (T - 1)]);
  int shock_term = 0;
  int generate_shock_free = 0;
  
  int n_below_threshold = 0;
  for (i in 1 : (C * (T - 1))) {
    n_below_threshold += (abs(diff[i]) < outlier_threshold) ? 1 : 0;
  }
  array[n_below_threshold] int indices_below_threshold;
  
  {
    int index = 1;
    for (i in 1 : (C * (T - 1))) {
      if (abs(diff[i]) < outlier_threshold) {
        indices_below_threshold[index] = i;
        index += 1;
      }
    }
  }
}
parameters {
  real<lower=0> epsilon_sigma;
  
  vector[C] raw_Delta1;
  array[hierarchical] real mu_Delta1;
  array[hierarchical] real sigma_Delta1;
  vector[C] raw_Delta2;
  array[hierarchical] real mu_Delta2;
  array[hierarchical] real sigma_Delta2;
  vector[C] raw_Delta3;
  array[hierarchical] real mu_Delta3;
  array[hierarchical] real sigma_Delta3;
  vector[C] raw_Delta4;
  array[hierarchical] real mu_Delta4;
  array[hierarchical] real sigma_Delta4;
  vector[C] raw_k;
  array[hierarchical] real mu_k;
  array[hierarchical] real sigma_k;
  vector[C] raw_z;
  array[hierarchical] real mu_z;
  array[hierarchical] real sigma_z;
}
transformed parameters {
  matrix[C, T - 1] shock = rep_matrix(0, C, T - 1);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  vector[C] first_transition = rep_vector(0, C);
  vector[C] final_transition = rep_vector(0, C);
  
  vector[C] Delta1;
  
  if (hierarchical == 1) {
    Delta1 = mu_Delta1[1] + exp(sigma_Delta1[1]) * raw_Delta1;
  } else {
    Delta1 = raw_Delta1;
  }
  if (Delta1_constrain == 1) {
    Delta1 = inv_logit(Delta1) * (Delta1_upper - Delta1_lower) + Delta1_lower;
  }
  vector[C] Delta2;
  
  if (hierarchical == 1) {
    Delta2 = mu_Delta2[1] + exp(sigma_Delta2[1]) * raw_Delta2;
  } else {
    Delta2 = raw_Delta2;
  }
  if (Delta2_constrain == 1) {
    Delta2 = inv_logit(Delta2) * (Delta2_upper - Delta2_lower) + Delta2_lower;
  }
  vector[C] Delta3;
  
  if (hierarchical == 1) {
    Delta3 = mu_Delta3[1] + exp(sigma_Delta3[1]) * raw_Delta3;
  } else {
    Delta3 = raw_Delta3;
  }
  if (Delta3_constrain == 1) {
    Delta3 = inv_logit(Delta3) * (Delta3_upper - Delta3_lower) + Delta3_lower;
  }
  vector[C] Delta4;
  
  if (hierarchical == 1) {
    Delta4 = mu_Delta4[1] + exp(sigma_Delta4[1]) * raw_Delta4;
  } else {
    Delta4 = raw_Delta4;
  }
  if (Delta4_constrain == 1) {
    Delta4 = inv_logit(Delta4) * (Delta4_upper - Delta4_lower) + Delta4_lower;
  }
  vector[C] k;
  
  if (hierarchical == 1) {
    k = mu_k[1] + exp(sigma_k[1]) * raw_k;
  } else {
    k = raw_k;
  }
  if (k_constrain == 1) {
    k = inv_logit(k) * (k_upper - k_lower) + k_lower;
  }
  vector[C] z;
  
  if (hierarchical == 1) {
    z = mu_z[1] + exp(sigma_z[1]) * raw_z;
  } else {
    z = raw_z;
  }
  if (z_constrain == 1) {
    z = inv_logit(z) * (z_upper - z_lower) + z_lower;
  }
  
  transition_function = to_matrix(
                                  rate_double_logistic(
                                    to_vector(y[ : , 1 : (T - 1)]),
                                    rep_vector_times(Delta1, T - 1),
                                    rep_vector_times(Delta2, T - 1),
                                    rep_vector_times(Delta3, T - 1),
                                    rep_vector_times(Delta4, T - 1),
                                    rep_vector_times(k, T - 1),
                                    rep_vector_times(z, T - 1)),
                                  C, T - 1);
  
  if (include_prior == 1) {
    first_transition = rate_double_logistic(rep_vector(grid[1], C), Delta1,
                         Delta2, Delta3, Delta4, k, z);
    final_transition = rate_double_logistic(rep_vector(grid[num_grid], C),
                         Delta1, Delta2, Delta3, Delta4, k, z);
  }
}
model {
  if (include_prior == 1) {
    to_vector(first_transition) ~ normal(0, 25);
    to_vector(final_transition) ~ normal(1.15 / 10, 0.5);
  }
  
  epsilon_sigma ~ std_normal();
  
  if (outlier_threshold < 1000) {
    diff[indices_below_threshold] ~ normal(
                                           to_vector(transition_function)[indices_below_threshold],
                                           epsilon_sigma);
  } else {
    diff ~ normal(to_vector(transition_function), epsilon_sigma);
  }
  
  if (hierarchical) {
    raw_Delta1 ~ std_normal();
    mu_Delta1[1] ~ normal(Delta1_prior_mean, Delta1_prior_sd);
    
    sigma_Delta1[1] ~ normal(-1, 0.5);
  } else {
    raw_Delta1 ~ normal(Delta1_prior_mean, Delta1_prior_sd);
  }
  if (hierarchical) {
    raw_Delta2 ~ std_normal();
    mu_Delta2[1] ~ normal(Delta2_prior_mean, Delta2_prior_sd);
    
    sigma_Delta2[1] ~ normal(-1, 0.5);
  } else {
    raw_Delta2 ~ normal(Delta2_prior_mean, Delta2_prior_sd);
  }
  if (hierarchical) {
    raw_Delta3 ~ std_normal();
    mu_Delta3[1] ~ normal(Delta3_prior_mean, Delta3_prior_sd);
    
    sigma_Delta3[1] ~ normal(-1, 0.5);
  } else {
    raw_Delta3 ~ normal(Delta3_prior_mean, Delta3_prior_sd);
  }
  if (hierarchical) {
    raw_Delta4 ~ std_normal();
    mu_Delta4[1] ~ normal(Delta4_prior_mean, Delta4_prior_sd);
    
    sigma_Delta4[1] ~ normal(-1, 0.5);
  } else {
    raw_Delta4 ~ normal(Delta4_prior_mean, Delta4_prior_sd);
  }
  if (hierarchical) {
    raw_k ~ std_normal();
    mu_k[1] ~ normal(k_prior_mean, k_prior_sd);
    
    sigma_k[1] ~ normal(-1, 0.5);
  } else {
    raw_k ~ normal(k_prior_mean, k_prior_sd);
  }
  if (hierarchical) {
    raw_z ~ std_normal();
    mu_z[1] ~ normal(z_prior_mean, z_prior_sd);
    
    sigma_z[1] ~ normal(-1, 0.5);
  } else {
    raw_z ~ normal(z_prior_mean, z_prior_sd);
  }
}
generated quantities {
  matrix[C, Tpred] eta;
  
  matrix[generate_shock_free * C, generate_shock_free * Tpred] eta_shockfree;
  matrix[shock_term * C, shock_term * (Tpred - 1)] shock2;
  if (shock_term == 1) 
    shock2 = rep_matrix(0, C, Tpred - 1);
  
  eta[1 : C, 1 : T] = y;
  if (shock_term == 1) {
    shock2[1 : C, 1 : (T - 1)] = shock;
  }
  
  if (generate_shock_free == 1) {
    eta_shockfree[1 : C, 1 : T] = y;
  }
  
  matrix[C, num_grid] transition_function_pred;
  vector[num_grid * hierarchical] transition_function_pred_mean;
  
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
  
  for (t in (T + 1) : Tpred) {
    vector[C] transition = rate_double_logistic(eta[ : , t - 1], Delta1,
                             Delta2, Delta3, Delta4, k, z);
    for (c in 1 : C) {
      eta[c, t] = eta[c, t - 1] + transition[c]
                  + error_rng(eta[c : c, t - 1], epsilon_params, 1);
      if (shock_term == 1) 
        eta[c, t] += shock2[c, t - 1];
    }
  }
  
  if (generate_shock_free == 1) {
    for (t in (T + 1) : Tpred) {
      vector[C] transition = rate_double_logistic(eta_shockfree[ : , 
                               t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
      for (c in 1 : C) {
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition[c]
                              + error_rng(eta_shockfree[c : c, t - 1],
                                          epsilon_params, 0);
      }
    }
  }
  
  for (i in 1 : num_grid) {
    transition_function_pred[ : , i] = rate_double_logistic(
                                         rep_vector(grid[i], C), Delta1,
                                         Delta2, Delta3, Delta4, k, z);
  }
  
  if (hierarchical == 1) {
    transition_function_pred_mean = rate_double_logistic(grid,
                                      rep_vector(
                                                 inv_logit(mu_Delta1[1])
                                                 * (Delta1_upper
                                                    - Delta1_lower)
                                                 + Delta1_lower, num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta2[1])
                                                 * (Delta2_upper
                                                    - Delta2_lower)
                                                 + Delta2_lower, num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta3[1])
                                                 * (Delta3_upper
                                                    - Delta3_lower)
                                                 + Delta3_lower, num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta4[1])
                                                 * (Delta4_upper
                                                    - Delta4_lower)
                                                 + Delta4_lower, num_grid),
                                      rep_vector(
                                                 inv_logit(mu_k[1])
                                                 * (k_upper - k_lower)
                                                 + k_lower, num_grid),
                                      rep_vector(
                                                 inv_logit(mu_z[1])
                                                 * (z_upper - z_lower)
                                                 + z_lower, num_grid));
  }
}

