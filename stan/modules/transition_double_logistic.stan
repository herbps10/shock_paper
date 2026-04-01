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
transformed parameters {
  transition_function = to_matrix(
    rate_double_logistic(
      to_vector(y[, 1:(T - 1)]),
      rep_vector_times(Delta[, 1], T - 1),
      rep_vector_times(Delta[, 2], T - 1),
      rep_vector_times(Delta[, 3], T - 1),
      rep_vector_times(Delta[, 4], T - 1),
      rep_vector_times(Delta[, 5], T - 1),
      rep_vector_times(Delta[, 6], T - 1)
    ), C, T - 1);
  
  if(include_prior == 1) {
    first_transition[1] = rate_double_logistic(rep_vector(grid[1], C), Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
    if(intermediate_grid_index > 0) {
      intermediate_transition[1] = rate_double_logistic(rep_vector(grid[intermediate_grid_index], C), Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
    }
    final_transition[1] = rate_double_logistic(rep_vector(grid[num_grid], C), Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
  }
}
generated quantities {
  // With shocks
  for(t in (T + 1):Tpred) {
    vector[C] transition = rate_double_logistic(eta[, t - 1], Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
    for(c in 1:C) {
      eta[c, t] = eta[c, t - 1] + transition[c] + error_rng(eta[c:c, t - 1], epsilon_params, 1);
      if(shock_term == 1) {
        if(shock_diff_mode == 1) {
          eta[c, t] += shock2[c, t] - shock2[c, t - 1];
        }
        else {
          eta[c, t] += shock2[c, t - 1];
        }
      }
    }
  }

  // Without shocks
  if(generate_shock_free == 1) {
    for(t in (T + 1):Tpred) {
      vector[C] transition = rate_double_logistic(eta_shockfree[, t - 1], Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
      for(c in 1:C) {
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition[c] + error_rng(eta_shockfree[c:c, t - 1], epsilon_params, 0);
      }
    }
  }

  for(i in 1:num_grid) {
    transition_function_pred[, i] = rate_double_logistic(rep_vector(grid[i], C), Delta[, 1], Delta[, 2], Delta[, 3], Delta[, 4], Delta[, 5], Delta[, 6]);
  }

  if(hierarchical == 1) {
    transition_function_pred_mean = rate_double_logistic(
      grid,
      rep_vector(inv_logit(mu_Delta[1][1]) * (Delta_upper[1] - Delta_lower[1]) + Delta_lower[1], num_grid), 
      rep_vector(inv_logit(mu_Delta[1][2]) * (Delta_upper[2] - Delta_lower[2]) + Delta_lower[2], num_grid),
      rep_vector(inv_logit(mu_Delta[1][3]) * (Delta_upper[3] - Delta_lower[3]) + Delta_lower[3], num_grid), 
      rep_vector(inv_logit(mu_Delta[1][4]) * (Delta_upper[4] - Delta_lower[4]) + Delta_lower[4], num_grid),
      rep_vector(inv_logit(mu_Delta[1][5]) * (Delta_upper[5] - Delta_lower[5]) + Delta_lower[5], num_grid),
      rep_vector(inv_logit(mu_Delta[1][6]) * (Delta_upper[6] - Delta_lower[6]) + Delta_lower[6], num_grid)
    );
  }
}
