functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
  
  real deboor(real x, vector ext_knots, row_vector a, int degree) {
    int k = degree + 1;
    row_vector[degree + 1] d;
    int n_ext_knots = rows(ext_knots);
    
    if (x <= ext_knots[1]) 
      return a[1];
    if (x >= ext_knots[n_ext_knots]) 
      return a[cols(a)];
    
    while (!(ext_knots[k + 1] > x) && k < n_ext_knots - degree - 1) {
      k = k + 1;
    }
    
    d = a[(k - degree) : (k)];
    
    for (r in 2 : (degree + 1)) {
      for (j in (k + r - degree - 1) : k) {
        int j2 = (k + r - degree - 1) - j + k;
        real alpha = (x - ext_knots[j2])
                     / (ext_knots[j2 + 1 + degree - (r - 1)] - ext_knots[j2]);
        
        d[j2 + 1 - k + degree] = (1 - alpha) * d[j2 - k + degree]
                                 + alpha * d[j2 + 1 - k + degree];
      }
    }
    
    return d[degree + 1];
  }
  
  real rate_spline(real P, real P_tilde, real P_tilde2, row_vector alpha,
                   vector ext_knots, int num_basis, int spline_degree) {
    return deboor((P - P_tilde) / (P_tilde2 - P_tilde), ext_knots, alpha,
                  spline_degree);
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
  
  int<lower=0, upper=1> alpha_constrain;
  real alpha_lower;
  real alpha_upper;
  real alpha_prior_mean;
  real alpha_prior_sd;
  
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  matrix[num_knots + spline_degree - 1, num_grid] B;
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
  print(intermediate_grid_index);
  
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
  
  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  
  ext_knots[1 : spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1) : (num_knots + 2 * spline_degree)] = rep_vector(
                                                                    knots[num_knots],
                                                                    spline_degree);
  ext_knots[(spline_degree + 1) : (num_knots + spline_degree)] = knots;
  
  real P_tilde = 5;
  real P_tilde2 = 110;
}
parameters {
  real<lower=0> epsilon_sigma;
  
  matrix[C, num_basis] raw_alpha;
  array[hierarchical] vector[num_basis] mu_alpha;
  array[hierarchical] vector<lower=0>[num_basis] sigma_alpha;
}
transformed parameters {
  matrix[C, T - 1] shock = rep_matrix(0, C, T - 1);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  array[include_prior] vector[C] first_transition;
  array[include_prior] vector[C] intermediate_transition;
  array[include_prior] vector[C] final_transition;
  
  if (include_prior == 1) {
    first_transition[1] = rep_vector(0, C);
    intermediate_transition[1] = rep_vector(0, C);
    final_transition[1] = rep_vector(0, C);
  }
  
  matrix[C, num_basis] alpha;
  
  for (n in 1 : num_basis) {
    if (hierarchical == 1) {
      alpha[ : , n] = mu_alpha[1][n] + sigma_alpha[1][n] * raw_alpha[ : , n];
    } else {
      alpha[ : , n] = raw_alpha[ : , n];
    }
    
    if (alpha_constrain == 1) {
      alpha[ : , n] = inv_logit(alpha[ : , n]) * (alpha_upper - alpha_lower)
                      + alpha_lower;
    }
  }
  
  for (c in 1 : C) {
    for (t in 1 : (T - 1)) {
      transition_function[c, t] = rate_spline(y[c, t], P_tilde, P_tilde2,
                                              alpha[c,  : ], ext_knots,
                                              num_basis, spline_degree);
    }
  }
  
  if (include_prior == 1) {
    for (c in 1 : C) 
      first_transition[1][c] = rate_spline(grid[1] / 110, 0, 1, alpha[c],
                                           ext_knots, num_basis,
                                           spline_degree);
    if (intermediate_grid_index > 0) {
      for (c in 1 : C) 
        intermediate_transition[1][c] = rate_spline(
                                                    grid[intermediate_grid_index]
                                                    / 110, 0, 1, alpha[c],
                                                    ext_knots, num_basis,
                                                    spline_degree);
    }
    for (c in 1 : C) 
      final_transition[1][c] = rate_spline(grid[num_grid] / 110, 0, 1,
                                           alpha[c], ext_knots, num_basis,
                                           spline_degree);
  }
}
model {
  if (include_prior == 1) {
    to_vector(first_transition[1]) ~ normal(0, 25);
    to_vector(intermediate_transition[1]) ~ normal(0, 4);
    to_vector(final_transition[1]) ~ normal(1.15 / 10, 0.5);
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
    to_vector(raw_alpha) ~ std_normal();
    mu_alpha[1] ~ normal(alpha_prior_mean, alpha_prior_sd);
    sigma_alpha[1] ~ std_normal();
  } else {
    to_vector(raw_alpha) ~ normal(alpha_prior_mean, alpha_prior_sd);
  }
  
  to_vector(alpha) ~ std_normal();
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
    for (c in 1 : C) {
      real transition = rate_spline(eta[c, t - 1], P_tilde, P_tilde2,
                                    alpha[c], ext_knots, num_basis,
                                    spline_degree);
      eta[c, t] = eta[c, t - 1] + transition
                  + error_rng(eta[c : c, t - 1], epsilon_params, 1);
      if (shock_term == 1) 
        eta[c, t] += shock2[c, t - 1];
    }
  }
  
  if (shock_term == 1) {
    for (t in (T + 1) : Tpred) {
      for (c in 1 : C) {
        real transition = rate_spline(eta_shockfree[c, t - 1], P_tilde,
                                      P_tilde2, alpha[c], ext_knots,
                                      num_basis, spline_degree);
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition
                              + error_rng(eta_shockfree[c : c, t - 1],
                                          epsilon_params, 0);
      }
    }
  }
  
  for (c in 1 : C) {
    for (i in 1 : num_grid) {
      transition_function_pred[c, i] = rate_spline(grid[i] / 110, 0, 1,
                                                   alpha[c], ext_knots,
                                                   num_basis, spline_degree);
    }
  }
  
  if (hierarchical == 1) {
    row_vector[num_basis] mu_alpha_pred = to_row_vector(mu_alpha[1]);
    if (alpha_constrain == 1) {
      mu_alpha_pred = inv_logit(mu_alpha_pred) * (alpha_upper - alpha_lower)
                      + alpha_lower;
    }
    for (i in 1 : num_grid) {
      transition_function_pred_mean[i] = rate_spline(grid[i] / 110, 0, 1,
                                                     mu_alpha_pred,
                                                     ext_knots, num_basis,
                                                     spline_degree);
    }
  }
}

