functions {
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
  real shock_rng(real nu_local, real c_slab, real tau,
                 int constrain_negative) {
    real shock_raw_pred;
    if (constrain_negative == 1) {
      shock_raw_pred = normal_lub_rng(0, 1, negative_infinity(), 0);
    } else {
      shock_raw_pred = normal_rng(0, 1);
    }
    real local_shrinkage_pred = student_t_rng(nu_local, 0, 1);
    real truncated_local_shrinkage_pred = sqrt(
                                               square(c_slab)
                                               * square(local_shrinkage_pred)
                                               ./ (square(c_slab)
                                                   + square(tau)
                                                     * square(
                                                              local_shrinkage_pred)));
    return shock_raw_pred * truncated_local_shrinkage_pred * tau;
  }
  
  real error_rng(vector x, vector epsilon_params, int shock) {
    real inner_sigma = epsilon_params[2];
    if (shock == 1) {
      real gamma = epsilon_params[1];
      real outer_sigma = epsilon_params[3];
      
      if (bernoulli_rng(gamma)) {
        return normal_rng(0, inner_sigma);
      } else {
        return normal_rng(0, outer_sigma);
      }
    } else {
      return normal_rng(0, inner_sigma);
    }
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
  
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  int<lower=0, upper=1> constrain_negative;
  
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
  
  shock_term = 1;
  generate_shock_free = 1;
  real tau = scale_global;
  int nu_local = 3;
  
  generate_shock_free = 1;
}
parameters {
  vector<upper=(constrain_negative == 1 ? 0 : positive_infinity())>[C * T] shock_raw;
  vector<lower=0>[C * T] lambda;
  real<lower=0> caux;
  
  real<lower=0> inner_epsilon_sigma;
  real<lower=inner_epsilon_sigma> outer_epsilon_sigma;
  real<lower=0, upper=1> gamma;
  
  matrix[C, D] raw_Delta;
  array[hierarchical] vector[D] mu_Delta;
  array[hierarchical] vector<lower=0>[D] sigma_Delta;
  array[hierarchical] cholesky_factor_corr[D] L_Omega_Delta;
}
transformed parameters {
  matrix[C, T] shock = rep_matrix(0, C, T);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  array[include_prior] vector[C] first_transition;
  array[include_prior] vector[C] intermediate_transition;
  array[include_prior] vector[C] final_transition;
  
  if (include_prior == 1) {
    first_transition[1] = rep_vector(0, C);
    intermediate_transition[1] = rep_vector(0, C);
    final_transition[1] = rep_vector(0, C);
  }
  
  real<lower=0> c_slab = slab_scale * sqrt(caux);
  vector<lower=0>[C * T] lambda_tilde = sqrt(
                                             c_slab ^ 2 * square(lambda)
                                             ./ (c_slab ^ 2
                                                 + tau ^ 2 * square(lambda)));
  shock = to_matrix(shock_raw .* lambda_tilde * tau * epsilon_sigma, C, T);
  
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
  
  shock_raw ~ std_normal();
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  lambda ~ student_t(nu_local, 0, 1);
  
  inner_epsilon_sigma ~ normal(0, 0.75);
  outer_epsilon_sigma ~ normal(0, 10);
  gamma ~ beta(1, 50);
  
  target += log_mix(gamma,
                    normal_lpdf(
                                to_vector(diff)
                                - to_vector(transition_function)
                                - to_vector(shock) | 0, inner_epsilon_sigma),
                    normal_lpdf(
                                to_vector(diff)
                                - to_vector(transition_function)
                                - to_vector(shock) | 0, outer_epsilon_sigma));
  
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
  matrix[shock_term * C, shock_term * Tpred] shock2;
  if (shock_term == 1) 
    shock2 = rep_matrix(0, C, Tpred);
  
  eta[1 : C, 1 : T] = y;
  if (shock_term == 1) {
    shock2[1 : C, 1 : T] = shock;
  }
  
  if (generate_shock_free == 1) {
    eta_shockfree[1 : C, 1 : T] = y - shock;
  }
  
  matrix[C, num_grid] transition_function_pred;
  vector[num_grid * hierarchical] transition_function_pred_mean;
  
  real lambda_tilde_sd = sd(lambda_tilde);
  
  for (t in T : Tpred) {
    for (c in 1 : C) {
      shock2[c, t] = shock_rng(nu_local, c_slab, tau, constrain_negative);
    }
  }
  
  matrix[C, (T - 1)] likelihood_ratio;
  for (c in 1 : C) {
    for (t in 1 : (T - 1)) {
      likelihood_ratio[c, t] = normal_lpdf(
                                           y[c, t + 1] - y[c, t]
                                           - transition_function[c, t]
                                           - shock[c, t] | 0,
                                           outer_epsilon_sigma)
                               - normal_lpdf(
                                             y[c, t + 1] - y[c, t]
                                             - transition_function[c, t]
                                             - shock[c, t] | 0,
                                             inner_epsilon_sigma);
    }
  }
  
  vector[3] epsilon_params;
  epsilon_params[1] = gamma;
  epsilon_params[2] = inner_epsilon_sigma;
  epsilon_params[3] = outer_epsilon_sigma;
  
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
      if (shock_term == 1) 
        eta[c, t] += shock2[c, t] - shock2[c, t - 1];
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

