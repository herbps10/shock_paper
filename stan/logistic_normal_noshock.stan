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
  int<lower=0, upper=1> shock_diff_mode;
  
  real<lower=0> epsilon_sigma_prior_mu;
  real<lower=0> epsilon_sigma_prior_sd;
  
  int<lower=1> D;
  
  array[D] int<lower=0, upper=1> Delta_constrain;
  vector[D] Delta_lower;
  vector[D] Delta_upper;
  vector[D] Delta_prior_mean;
  vector[D] Delta_prior_sd;
}
transformed data {
  vector[C * (T - 1)] diff = to_vector(y[ : , 2 : T] - y[ : , 1 : (T - 1)]);
  int shock_term = 0;
  int generate_shock_free = 0;
  
  int intermediate_grid_index = 0;
  for (i in 1 : num_grid) {
    if (grid[i] == 90) {
      intermediate_grid_index = i;
      break;
    }
  }
  
  int T_shocks = T;
  if (shock_diff_mode == 0) {
    T_shocks = T - 1;
  }
}
parameters {
  real<lower=0> epsilon_sigma;
  
  matrix[C, D] raw_Delta;
  array[hierarchical] vector[D] mu_Delta;
  array[hierarchical] vector<lower=0>[D] sigma_Delta;
  array[hierarchical] cholesky_factor_corr[D] L_Omega_Delta;
}
transformed parameters {
  matrix[C, T_shocks] shock = rep_matrix(0, C, T_shocks);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  array[include_prior] vector[C] first_transition;
  array[include_prior] vector[C] intermediate_transition;
  array[include_prior] vector[C] final_transition;
  
  if (include_prior == 1) {
    first_transition[1] = rep_vector(0, C);
    intermediate_transition[1] = rep_vector(0, C);
    final_transition[1] = rep_vector(0, C);
  }
  
  matrix[C, D] Delta;
  
  if (hierarchical == 1) {
    Delta = rep_matrix(mu_Delta[1]', C)
            + (diag_pre_multiply(sigma_Delta[1], L_Omega_Delta[1])
               * raw_Delta')';
  } else {
    Delta = raw_Delta;
  }
  
  for (n in 1 : D) {
    if (Delta_constrain[n] == 1) {
      Delta[ : , n] = inv_logit(Delta[ : , n])
                      * (Delta_upper[n] - Delta_lower[n]) + Delta_lower[n];
    }
  }
  
  transition_function = to_matrix(
                                  rate_double_logistic(
                                    to_vector(y[ : , 1 : (T - 1)]),
                                    rep_vector_times(Delta[ : , 1], T - 1),
                                    rep_vector_times(Delta[ : , 2], T - 1),
                                    rep_vector_times(Delta[ : , 3], T - 1),
                                    rep_vector_times(Delta[ : , 4], T - 1),
                                    rep_vector_times(Delta[ : , 5], T - 1),
                                    rep_vector_times(Delta[ : , 6], T - 1)),
                                  C, T - 1);
  
  if (include_prior == 1) {
    first_transition[1] = rate_double_logistic(rep_vector(grid[1], C),
                            Delta[ : , 1], Delta[ : , 2], Delta[ : , 3],
                            Delta[ : , 4], Delta[ : , 5], Delta[ : , 6]);
    if (intermediate_grid_index > 0) {
      intermediate_transition[1] = rate_double_logistic(
                                     rep_vector(
                                                grid[intermediate_grid_index],
                                                C),
                                     Delta[ : , 1], Delta[ : , 2],
                                     Delta[ : , 3], Delta[ : , 4],
                                     Delta[ : , 5], Delta[ : , 6]);
    }
    final_transition[1] = rate_double_logistic(rep_vector(grid[num_grid], C),
                            Delta[ : , 1], Delta[ : , 2], Delta[ : , 3],
                            Delta[ : , 4], Delta[ : , 5], Delta[ : , 6]);
  }
}
model {
  if (include_prior == 1) {
    to_vector(first_transition[1]) ~ normal(0, 25);
    to_vector(intermediate_transition[1]) ~ normal(0, 4);
    to_vector(final_transition[1]) ~ normal(1.15 / 10, 0.5);
  }
  
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);
  if (shock_diff_mode == 1) {
    diff ~ normal(
                  to_vector(
                            transition_function + shock[ : , 2 : T]
                            - shock[ : , 1 : (T - 1)]),
                  epsilon_sigma);
  } else {
    diff ~ normal(to_vector(transition_function + shock), epsilon_sigma);
  }
  
  if (hierarchical) {
    to_vector(raw_Delta) ~ std_normal();
    mu_Delta[1] ~ normal(Delta_prior_mean, Delta_prior_sd);
    sigma_Delta[1] ~ std_normal();
    L_Omega_Delta[1] ~ lkj_corr_cholesky(1.0);
  } else {
    for (d in 1 : D) {
      to_vector(raw_Delta[ : , d]) ~ normal(Delta_prior_mean[d],
                                            Delta_prior_sd[d]);
    }
  }
}
generated quantities {
  matrix[C, Tpred] eta;
  
  matrix[generate_shock_free * C, generate_shock_free * Tpred] eta_shockfree;
  matrix[shock_term * C, shock_term * (T_shocks + Tpred - T)] shock2;
  if (shock_term == 1) 
    shock2 = rep_matrix(0, C, T_shocks + Tpred - T);
  
  eta[1 : C, 1 : T] = y;
  if (shock_term == 1) {
    shock2[1 : C, 1 : T_shocks] = shock;
  }
  
  if (generate_shock_free == 1) {
    if (shock_diff_mode == 0) {
      eta_shockfree[ : , 1] = rep_vector(0, C);
      eta_shockfree[ : , 2 : T] = y[ : , 2 : T] - shock;
    } else {
      eta_shockfree[ : , 1 : T] = y - shock;
    }
  }
  
  matrix[C, num_grid] transition_function_pred;
  vector[num_grid * hierarchical] transition_function_pred_mean;
  
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
  
  corr_matrix[D * hierarchical] Omega_Delta;
  if (hierarchical) 
    Omega_Delta = multiply_lower_tri_self_transpose(L_Omega_Delta[1]);
  
  for (t in (T + 1) : Tpred) {
    vector[C] transition = rate_double_logistic(eta[ : , t - 1],
                             Delta[ : , 1], Delta[ : , 2], Delta[ : , 3],
                             Delta[ : , 4], Delta[ : , 5], Delta[ : , 6]);
    for (c in 1 : C) {
      eta[c, t] = eta[c, t - 1] + transition[c]
                  + error_rng(eta[c : c, t - 1], epsilon_params, 1);
      if (shock_term == 1) {
        if (shock_diff_mode == 1) {
          eta[c, t] += shock2[c, t] - shock2[c, t - 1];
        } else {
          eta[c, t] += shock2[c, t - 1];
        }
      }
    }
  }
  
  if (generate_shock_free == 1) {
    for (t in (T + 1) : Tpred) {
      vector[C] transition = rate_double_logistic(eta_shockfree[ : , 
                               t - 1], Delta[ : , 1], Delta[ : , 2],
                               Delta[ : , 3], Delta[ : , 4], Delta[ : , 5],
                               Delta[ : , 6]);
      for (c in 1 : C) {
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition[c]
                              + error_rng(eta_shockfree[c : c, t - 1],
                                          epsilon_params, 0);
      }
    }
  }
  
  for (i in 1 : num_grid) {
    transition_function_pred[ : , i] = rate_double_logistic(
                                         rep_vector(grid[i], C),
                                         Delta[ : , 1], Delta[ : , 2],
                                         Delta[ : , 3], Delta[ : , 4],
                                         Delta[ : , 5], Delta[ : , 6]);
  }
  
  if (hierarchical == 1) {
    transition_function_pred_mean = rate_double_logistic(grid,
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][1])
                                                 * (Delta_upper[1]
                                                    - Delta_lower[1])
                                                 + Delta_lower[1], num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][2])
                                                 * (Delta_upper[2]
                                                    - Delta_lower[2])
                                                 + Delta_lower[2], num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][3])
                                                 * (Delta_upper[3]
                                                    - Delta_lower[3])
                                                 + Delta_lower[3], num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][4])
                                                 * (Delta_upper[4]
                                                    - Delta_lower[4])
                                                 + Delta_lower[4], num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][5])
                                                 * (Delta_upper[5]
                                                    - Delta_lower[5])
                                                 + Delta_lower[5], num_grid),
                                      rep_vector(
                                                 inv_logit(mu_Delta[1][6])
                                                 * (Delta_upper[6]
                                                    - Delta_lower[6])
                                                 + Delta_lower[6], num_grid));
  }
}

