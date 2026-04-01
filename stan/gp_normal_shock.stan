functions {
  real flambda(real L, int m) {
    real lam;
    lam = ((m * pi()) / (2 * L)) ^ 2;
    
    return lam;
  }
  real spd(real alpha, real rho, real w) {
    return (alpha ^ 2) * sqrt(2 * pi()) * rho
           * exp(-0.5 * (rho ^ 2) .* (w ^ 2));
  }
  
  vector phi(real L, int m, vector x) {
    vector[rows(x)] fi;
    fi = 1 / sqrt(L) * sin(m * pi() / (2 * L) * (x + L));
    
    return fi;
  }
  
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
  real shock_rng(real nu_local, real c_slab, real tau) {
    real shock_raw_pred = normal_lub_rng(0, 1, negative_infinity(), 0);
    
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
    return normal_rng(0, epsilon_params[1]);
  }
  
  vector constraint(vector x, vector ys) {
    return inv_logit(x) * 30;
  }
  
  vector rate_gp(vector x, matrix PHI, vector beta) {
    return constraint(PHI * beta, x);
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
  
  real L;
  int<lower=1> M;
  
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  
  real<lower=0> epsilon_sigma_prior_mu;
  real<lower=0> epsilon_sigma_prior_sd;
  
  array[M] int<lower=0, upper=1> beta_constrain;
  vector[M] beta_lower;
  vector[M] beta_upper;
  vector[M] beta_prior_mean;
  vector[M] beta_prior_sd;
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
  
  matrix[C * (T - 1), M] PHI;
  matrix[num_grid, M] PHI_grid;
  for (m in 1 : M) {
    PHI_grid[ : , m] = phi(L, m, grid / 110);
    for (c in 1 : C) {
      int start = (c - 1) * (T - 1) + 1;
      int end = c * (T - 1);
      PHI[start : end, m] = phi(L, m, to_vector((y[c, 1 : (T - 1)]) / 110));
    }
  }
  
  shock_term = 1;
  generate_shock_free = 1;
  real tau = scale_global;
  int nu_local = 3;
}
parameters {
  vector<upper=0>[C * T] shock_raw;
  vector<lower=0>[C * T] lambda;
  real<lower=0> caux;
  
  real<lower=0> epsilon_sigma;
  
  matrix[C, M] raw_beta;
  array[hierarchical] vector[M] mu_beta;
  array[hierarchical] vector<lower=0>[M] sigma_beta;
  
  real<lower=0> rho;
  real<lower=0> alpha;
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
  
  matrix[C, M] beta;
  
  for (n in 1 : M) {
    if (hierarchical == 1) {
      beta[ : , n] = mu_beta[1][n] + sigma_beta[1][n] * raw_beta[ : , n];
    } else {
      beta[ : , n] = raw_beta;
    }
    
    if (beta_constrain[n] == 1) {
      beta[ : , n] = inv_logit(beta[ : , n])
                     * (beta_upper[n] - beta_lower[n]) + beta_lower[n];
    }
  }
  
  row_vector[M] diagSPD;
  matrix[C, M] SPD_beta;
  
  for (m in 1 : M) {
    diagSPD[m] = sqrt(spd(alpha, rho, sqrt(flambda(L, m))));
  }
  
  for (c in 1 : C) {
    SPD_beta[c,  : ] = diagSPD .* beta[c,  : ];
    int start = (c - 1) * (T - 1) + 1;
    int end = c * (T - 1);
    transition_function[c,  : ] = to_row_vector(
                                                rate_gp(
                                                        to_vector(
                                                                  y[c, 1 : (
                                                                  T - 1)]),
                                                        PHI[start : end,  : ],
                                                        to_vector(
                                                                  SPD_beta[c,  : ])));
  }
  
  if (include_prior == 1) {
    for (c in 1 : C) 
      first_transition[1][c] = rate_gp(grid[1 : 1],
                                       to_matrix(PHI_grid[1,  : ], 1, M),
                                       to_vector(SPD_beta[c,  : ]))[1];
    if (intermediate_grid_index > 0) {
      for (c in 1 : C) 
        intermediate_transition[1][c] = rate_gp(
                                                grid[intermediate_grid_index : intermediate_grid_index],
                                                to_matrix(PHI_grid[1,  : ],
                                                          1, M),
                                                to_vector(SPD_beta[c,  : ]))[1];
    }
    for (c in 1 : C) 
      final_transition[1][c] = rate_gp(grid[num_grid : num_grid],
                                       to_matrix(PHI_grid[num_grid,  : ], 1,
                                                 M),
                                       to_vector(SPD_beta[c,  : ]))[1];
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
  
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);
  to_vector(y[ : , 2 : T] - y[ : , 1 : (T - 1)]) ~ normal(
                                                          to_vector(
                                                                    transition_function
                                                                    + shock[ : , 2 : T]
                                                                    - shock[ : , 1 : (
                                                                    T - 1)]),
                                                          epsilon_sigma);
  
  if (hierarchical) {
    to_vector(raw_beta) ~ std_normal();
    mu_beta[1] ~ normal(beta_prior_mean, beta_prior_sd);
    sigma_beta[1] ~ std_normal();
  } else {
    to_vector(raw_beta) ~ normal(beta_prior_mean, beta_prior_sd);
  }
  
  alpha ~ std_normal();
  rho ~ inv_gamma(5, 5);
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
      shock2[c, t] = shock_rng(nu_local, c_slab, tau);
    }
  }
  
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
  
  for (t in (T + 1) : Tpred) {
    for (c in 1 : C) {
      matrix[1, M] PHI_pred;
      for (m in 1 : M) 
        PHI_pred[ : , m] = phi(L, m, eta[c : c, t - 1] / 110);
      real transition = rate_gp(eta[c : c, t - 1], PHI_pred,
                                to_vector(SPD_beta[c,  : ]))[1];
      
      eta[c, t] = eta[c, t - 1] + transition
                  + error_rng(eta[c : c, t - 1], epsilon_params, 1);
      
      if (shock_term == 1) 
        eta[c, t] += shock2[c, t - 1];
    }
  }
  
  if (generate_shock_free == 1) {
    for (t in (T + 1) : Tpred) {
      for (c in 1 : C) {
        matrix[1, M] PHI_pred;
        for (m in 1 : M) 
          PHI_pred[ : , m] = phi(L, m, eta_shockfree[c : c, t - 1] / 110);
        real transition = rate_gp(eta_shockfree[c : c, t - 1], PHI_pred,
                                  to_vector(SPD_beta[c,  : ]))[1];
        
        eta_shockfree[c, t] = eta_shockfree[c, t - 1] + transition
                              + error_rng(eta_shockfree[c : c, t - 1],
                                          epsilon_params, 0);
      }
    }
  }
  
  for (c in 1 : C) {
    transition_function_pred[c,  : ] = to_row_vector(
                                                     rate_gp(grid / 110,
                                                             PHI_grid,
                                                             to_vector
                                                             (
                                                              SPD_beta[c,  : ])));
  }
  
  if (hierarchical == 1) {
    vector[M] SPD_beta_mu = to_vector(diagSPD) .* mu_beta[1];
    transition_function_pred_mean = rate_gp(grid / 110, PHI_grid,
                                            SPD_beta_mu);
  }
}

