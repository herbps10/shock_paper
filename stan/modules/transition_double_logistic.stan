functions {
  vector rate_double_logistic(vector x, vector Delta1, vector Delta2, vector Delta3, vector Delta4, vector k, vector z) {
      real A1 = 4.4;
      real A2 = 0.5;
      return k .* inv(1 + exp(-A1 .* inv(Delta2) .* (x - Delta1 - A2 * Delta2))) + (z - k) .* inv(1 + exp(-A1 * inv(Delta4) .* (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
  }

  vector rep_vector_times(vector x, int times) {
    int s = size(x);
    vector[s  * times] y;
    for(i in 1:times) {
      y[((i - 1) * s + 1):(i * s)] = x;
    }
    return y;
  }
}
parameters {
}
transformed parameters {
  transition_function = to_matrix(
    rate_double_logistic(
      to_vector(y[, 1:(T - 1)]),
      rep_vector_times(Delta1, T - 1),
      rep_vector_times(Delta2, T - 1),
      rep_vector_times(Delta3, T - 1),
      rep_vector_times(Delta4, T - 1),
      rep_vector_times(k, T - 1),
      rep_vector_times(z, T - 1)
    ), C, T - 1);
  
  if(include_prior == 1) {
    first_transition[1] = rate_double_logistic(rep_vector(grid[1], C), Delta1, Delta2, Delta3, Delta4, k, z);
    if(intermediate_grid_index > 0) {
      intermediate_transition[1] = rate_double_logistic(rep_vector(grid[intermediate_grid_index], C), Delta1, Delta2, Delta3, Delta4, k, z);
    }
    final_transition[1] = rate_double_logistic(rep_vector(grid[num_grid], C), Delta1, Delta2, Delta3, Delta4, k, z);
  }
}
generated quantities {
  // With shocks
  for(t in (T + 1):Tpred) {
    vector[C] transition = rate_double_logistic(eta[, t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
    for(c in 1:C) {
      eta[c, t] = eta[c, t - 1] + transition[c] + error_rng(eta[c:c, t - 1], epsilon_params, 1);
      if(shock_term == 1) eta[c, t] += shock2[c, t - 1];
    }
  }

  // Without shocks
  if(generate_shock_free == 1) {
    for(t in (T + 1):Tpred) {
      vector[C] transition = rate_double_logistic(eta_shockfree[, t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
      for(c in 1:C) {
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition[c] + error_rng(eta_shockfree[c:c, t - 1], epsilon_params, 0);
      }
    }
  }

  for(i in 1:num_grid) {
    transition_function_pred[, i] = rate_double_logistic(rep_vector(grid[i], C), Delta1, Delta2, Delta3, Delta4, k, z);
  }

  if(hierarchical == 1) {
    transition_function_pred_mean = rate_double_logistic(
      grid,
      rep_vector(inv_logit(mu_Delta1[1]) * (Delta1_upper - Delta1_lower) + Delta1_lower, num_grid), 
      rep_vector(inv_logit(mu_Delta2[1]) * (Delta2_upper - Delta2_lower) + Delta2_lower, num_grid),
      rep_vector(inv_logit(mu_Delta3[1]) * (Delta3_upper - Delta3_lower) + Delta3_lower, num_grid), 
      rep_vector(inv_logit(mu_Delta4[1]) * (Delta4_upper - Delta4_lower) + Delta4_lower, num_grid),
      rep_vector(inv_logit(mu_k[1]) * (k_upper - k_lower) + k_lower, num_grid),
      rep_vector(inv_logit(mu_z[1]) * (z_upper - z_lower) + z_lower, num_grid)
    );
  }
}
