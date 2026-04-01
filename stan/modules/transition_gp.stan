functions {
  vector constraint(vector x, vector ys) {
	  return inv_logit(x) * 30;
	}

  vector rate_gp(vector x, matrix PHI, vector beta) {
    return constraint(PHI * beta, x);
  }
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
}
transformed parameters {
  row_vector[M] diagSPD;
  matrix[C, M] SPD_beta;

  for(m in 1:M) {
    diagSPD[m] = sqrt(spd(alpha, rho, sqrt(flambda(L, m))));
  }

  for(c in 1:C) {
    SPD_beta[c, ] = diagSPD .* beta[c, ];
    int start = (c - 1) * (T - 1) + 1;
    int end = c * (T - 1);
    transition_function[c, ] = to_row_vector(rate_gp(to_vector(y[c, 1:(T - 1)]), PHI[start:end, ], to_vector(SPD_beta[c, ])));
  }
  
  if(include_prior == 1) {
    for(c in 1:C) first_transition[1][c] = rate_gp(grid[1:1], to_matrix(PHI_grid[1, ], 1, M), to_vector(SPD_beta[c, ]))[1];
    if(intermediate_grid_index > 0) {
      for(c in 1:C) intermediate_transition[1][c] = rate_gp(grid[intermediate_grid_index:intermediate_grid_index], to_matrix(PHI_grid[1, ], 1, M), to_vector(SPD_beta[c, ]))[1];
    }
    for(c in 1:C) final_transition[1][c] = rate_gp(grid[num_grid:num_grid], to_matrix(PHI_grid[num_grid, ], 1, M), to_vector(SPD_beta[c, ]))[1];
  }
}
model {
  alpha ~ std_normal();
  rho ~ inv_gamma(5, 5);
}
generated quantities {
  // With shocks
  for(t in (T + 1):Tpred) {
    for(c in 1:C) {
      matrix[1, M] PHI_pred;
      for(m in 1:M) PHI_pred[, m] = phi(L, m, eta[c:c, t - 1] / 110);
      real transition = rate_gp(eta[c:c, t - 1], PHI_pred, to_vector(SPD_beta[c, ]))[1];

      eta[c, t] = eta[c, t - 1] + transition + error_rng(eta[c:c, t - 1], epsilon_params, 1);

      if(shock_term == 1) eta[c, t] += shock2[c, t - 1];
    }
  }

  // Without shocks
  if(generate_shock_free == 1) {
    for(t in (T + 1):Tpred) {
      for(c in 1:C) {
        matrix[1, M] PHI_pred;
        for(m in 1:M) PHI_pred[, m] = phi(L, m, eta_shockfree[c:c, t - 1] / 110);
        real transition = rate_gp(eta_shockfree[c:c, t - 1], PHI_pred, to_vector(SPD_beta[c, ]))[1];

        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition + error_rng(eta_shockfree[c:c, t - 1], epsilon_params, 0);
      }
    }
  }

  for(c in 1:C) {
    transition_function_pred[c, ] = to_row_vector(rate_gp(grid / 110, PHI_grid, to_vector(SPD_beta[c, ])));
  }

  if(hierarchical == 1) {
    vector[M] SPD_beta_mu = to_vector(diagSPD) .* mu_beta[1];
    transition_function_pred_mean = rate_gp(grid / 110, PHI_grid, SPD_beta_mu);
  }
}
