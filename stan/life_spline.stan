functions {
  #include ./scale_blocks.stan
  #include ./deboor.stan

  real rate_spline(real P, real P_tilde, real P_tilde2, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
    return deboor((P - P_tilde) / (P_tilde2 - P_tilde), ext_knots, a, spline_degree);
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int C; // Number of countries
  int R;
  int t_last;

  vector[N] y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=1> held_out;

  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  matrix[num_knots + spline_degree - 1, num_grid] B;

  int num_grid;
  vector[num_grid] grid;

  real a_lower_bound;
  real a_upper_bound;
  
}
transformed data {
  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;

  array[C] int final_observed = rep_array(0, C);

  matrix[C, t_last] ymat = rep_matrix(0, C, t_last);

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;

  for(i in 1:N) {
    ymat[country[i], time[i]] = y[i];

    if(held_out[i] == 0 && time[i] > final_observed[country[i]]) {
      final_observed[country[i]] = time[i];
    }
  }

  real P_tilde = 5;
  real P_tilde2 = 110;
  
  int hierarchical = 1;
}

parameters {
  // Spline rate vs. level function
  vector[hierarchical * (num_basis)] a_mu;
  vector<lower=0>[hierarchical * (num_basis)] a_sigma;
  matrix[C, num_basis] a_raw;

  real<lower=0> epsilon_variance;
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, num_basis] a = rep_matrix(0, C, num_basis);
  
  real epsilon_sd = sqrt(epsilon_variance);

  // Initialize the spline coefficients
  for(i in 1:(num_basis - 3)) {
    if(hierarchical == 1) {
      a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_mu[i] + a_raw[,i] * a_sigma[i]);
    }
    else {
      a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_raw[,i]);
    }
  }
  if(hierarchical == 1) {
    a[, num_basis - 2] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis - 2] + a_raw[,num_basis - 2] * a_sigma[num_basis - 2]);
    a[, num_basis - 1] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis - 1] + a_raw[,num_basis - 1] * a_sigma[num_basis - 1]);
    a[, num_basis]     = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis]     + a_raw[,num_basis] * a_sigma[num_basis]);
  }
  else {
    a[, num_basis - 2] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_raw[,num_basis - 2]);
    a[, num_basis - 1] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_raw[,num_basis - 1]);
    a[, num_basis]     = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_raw[,num_basis]);
  }

  for(c in 1:C) {
    for(t in 2:t_last) {
      transition_function[c, t] = rate_spline(ymat[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
    }
  }


}

model {
  if(hierarchical == 1) {
    a_mu ~ normal(0, 15);
    a_sigma ~ normal(0, 5); // increasing prior variance based on checks
  }
  to_vector(a_raw) ~ std_normal();

  epsilon_variance ~ inv_gamma(0.1, 0.1);


  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > 1) {
      (ymat[country[i], time[i]] - ymat[country[i], time[i] - 1]) ~ normal(transition_function[country[i], time[i]], epsilon_sd);
    }
  }
}
generated quantities {
  matrix[C, T] eta;
  matrix[C, num_grid] transition_function_pred;
  vector[num_grid] transition_function_mean;

  if(hierarchical == 1) {
    vector[num_basis] a_mean;
    a_mean[1:(num_basis - 3)] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_mu[1:(num_basis - 3)]);
    a_mean[num_basis - 2] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis - 2]);
    a_mean[num_basis - 1] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis - 1]);
    a_mean[num_basis] = a_lower_bound + (1.15/5.0 - a_lower_bound) * inv_logit(a_mu[num_basis]);

    for(i in 1:num_grid) {
      transition_function_mean[i] = rate_spline(grid[i], 0, 1, to_row_vector(a_mean), ext_knots, num_basis, spline_degree);
    }
  }

  for(c in 1:C) {
    eta[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];

    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_sd);
      real transition = rate_spline(eta[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
      eta[c, t] = eta[c, t - 1] + transition + error;
    }

    for(i in 1:num_grid) {
      transition_function_pred[c, i] = rate_spline(grid[i], 0, 1, a[c,], ext_knots, num_basis, spline_degree);
    }
  }
}
