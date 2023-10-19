functions {
  #include ./scale_blocks.stan
  #include ./deboor.stan

  real rate_spline(real P, real P_tilde, real P_tilde2, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
    return deboor((P - P_tilde) / P_tilde2, ext_knots, a, spline_degree);
  }

  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }

  real shock_rng(real nu_local, real c_slab, real global_shrinkage) {
    real shock_raw_pred = normal_lub_rng(0, 1, negative_infinity(), 0);
    //real shock_raw_pred = normal_rng(0, 1);
    real local_shrinkage_pred = student_t_rng(nu_local, 0, 1);
    real truncated_local_shrinkage_pred = sqrt(c_slab^2 * square(local_shrinkage_pred) ./ (c_slab^2 + global_shrinkage^2 * square(local_shrinkage_pred)));
    return shock_raw_pred * truncated_local_shrinkage_pred * global_shrinkage;
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int C; // Number of countries
  int t_last;

  array[N] y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=1> held_out;

  int num_knots;
  vector[num_knots] knots;

  int spline_degree;

  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;

  real a_lower_bound;
  real a_upper_bound;

  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
}
transformed data {
  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  int num_constrained_zero = spline_degree + 1;

  matrix[C, t_last] ymat = rep_matrix(0, C, t_last);

  array[C] int final_observed = rep_array(0, C);
  int n_shocks = 0;

  real<lower=1> nu_global = 1;
  real<lower=1> nu_local = 1;

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;

  for(i in 1:N) {
    ymat[country[i], time[i]] = y[i];

    if(held_out[i] == 0 && time[i] > final_observed[country[i]]) {
      final_observed[country[i]] = time[i];
    }
  }

  for(c in 1:C) {
    //n_shocks += final_observed[c] - 2;
    n_shocks += final_observed[c] - 1;
  }

  real P_tilde = 15;
  real P_tilde2 = 85;
}

parameters {
  // Spline rate vs. level function
  vector[num_basis - 2] a_mu;
  matrix[C, num_basis - 2] a_raw;
  vector<lower=0>[num_basis - 2] a_sigma;

  //vector[C] Omega_raw;
  //real P_tilde2_mu;
  //real<lower=0> P_tilde2_sigma;
  //vector[C] P_tilde2_raw;

  real<lower=0> epsilon_scale;

  vector<upper=0>[n_shocks] shock_raw;
  //vector[n_shocks] shock_raw;
  real<lower=0> global_shrinkage;
  vector<lower=0>[n_shocks] local_shrinkage; // called lambda in paper
  real<lower=0> caux;
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, t_last] gamma = rep_matrix(0, C, t_last);
  matrix[C, num_basis] a;

  matrix[C, t_last] shock = rep_matrix(0, C, t_last);

  //vector[C] P_tilde = rep_vector(15, C); // Lower asymptote
  //vector[C] P_tilde2 = 50 + inv_logit(P_tilde2_mu + P_tilde2_raw * P_tilde2_sigma) * 50; // Upper asymptote

  real<lower=0> c_slab = slab_scale * sqrt(caux);
  vector<lower=0>[n_shocks] truncated_local_shrinkage; // called lambda_tilde in paper

  {
    vector[n_shocks] shock_shrinkage;
    truncated_local_shrinkage = sqrt(c_slab^2 * square(local_shrinkage) ./ (c_slab^2 + global_shrinkage^2 * square(local_shrinkage)));
    shock_shrinkage = shock_raw .* truncated_local_shrinkage * global_shrinkage;

    int index = 1;
    for(c in 1:C) {
      int final_index = index + final_observed[c] - 2;
      shock[c, 1:(final_observed[c] - 1)] = to_row_vector(shock_shrinkage[index:final_index]);
      index += (final_observed[c] - 1);
    }
  }

  // Initialize the non-zero spline coefficients
  for(i in 1:(num_basis - 3)) {
    a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_mu[i] + a_raw[,i] * a_sigma[i]);
    //a[, i] = a_lower_bound + exp(a_mu[i] + a_raw[,i] * a_sigma[i]);
  }
  a[, num_basis - 2] = a_lower_bound + (1.15 - a_lower_bound) * inv_logit(a_mu[num_basis - 2] + a_raw[,num_basis - 2] * a_sigma[num_basis - 2]);

  for(c in 1:C) {
    a[c, (num_basis - 1):num_basis] = rep_row_vector(a[c, num_basis - 2], 2);
    for(t in 2:final_observed[c]) {
      transition_function[c, t] = rate_spline(ymat[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
      gamma[c, t] = transition_function[c, t] + shock[c, t] - shock[c, t - 1];
    }
  }
}

model {
  a_mu ~ normal(0, 3);
  to_vector(a_raw) ~ std_normal();
  a_sigma ~ std_normal();

  epsilon_scale ~ inv_gamma(0.1, 0.1);

  //P_tilde2_mu ~ std_normal();
  //P_tilde2_sigma ~ std_normal();
  //P_tilde2_raw ~ std_normal();

  shock_raw ~ std_normal();
  local_shrinkage ~ student_t(nu_local, 0, 1);
  global_shrinkage ~ student_t(nu_global, 0, scale_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > 1) {
      (ymat[country[i], time[i]] - ymat[country[i], time[i] - 1]) ~ normal(gamma[country[i], time[i]], epsilon_scale);
    }
  }
}
generated quantities {
  matrix[C, T] eta;
  matrix[C, T] eta_crisisfree;
  matrix[C, T] shock2;

  matrix[C, num_grid] transition_function_pred;
  vector[num_grid] transition_function_mean;

 {
    vector[num_basis] a_mean;
    a_mean[1:(num_basis - 3)] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_mu[1:(num_basis - 3)]);
    a_mean[num_basis - 2] = a_lower_bound + (1.15 - a_lower_bound) * inv_logit(a_mu[num_basis - 2]);
    a_mean[(num_basis - 1):num_basis] = rep_vector(a_mean[num_basis - 2], 2);
    for(i in 1:num_grid) {
      transition_function_mean[i] = rate_spline(grid[i], 0, 1, to_row_vector(a_mean), ext_knots, num_basis, spline_degree);
    }
  }

  for(c in 1:C) {
    eta_crisisfree[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];

    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition_crisisfree = rate_spline(eta_crisisfree[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
      eta_crisisfree[c, t] = eta_crisisfree[c, t - 1] + transition_crisisfree + error;
    }

    eta[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];
    shock2[c, 1:(final_observed[c])] = shock[c, 1:(final_observed[c])];

    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);

      shock2[c, t] = shock_rng(nu_local, c_slab, global_shrinkage);

      real transition = rate_spline(eta[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);

      eta[c, t] = eta[c, t - 1] + transition + error + shock2[c, t] - shock2[c, t - 1];
    }

    for(i in 1:num_grid) {
      transition_function_pred[c, i] = rate_spline(grid[i], 0, 1, a[c,], ext_knots, num_basis, spline_degree);
    }
  }
}
