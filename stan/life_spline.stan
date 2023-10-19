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
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, num_basis] a = rep_matrix(0, C, num_basis);

  // Initialize the spline coefficients
  for(i in 1:(num_basis - 3)) {
    a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(a_mu[i] + a_raw[,i] * a_sigma[i]);
  }
  a[, num_basis - 2] = a_lower_bound + (1.15 - a_lower_bound) * inv_logit(a_mu[num_basis - 2] + a_raw[,num_basis - 2] * a_sigma[num_basis - 2]);

  for(c in 1:C) {
    a[c, (num_basis - 1):num_basis] = rep_row_vector(a[c, num_basis - 2], 2);
    for(t in 2:t_last) {
      transition_function[c, t] = rate_spline(ymat[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
    }
  }
}

model {
  a_mu ~ std_normal();
  a_sigma ~ std_normal();
  to_vector(a_raw) ~ std_normal();

  epsilon_scale ~ inv_gamma(0.1, 0.1);

  //P_tilde2_mu ~ std_normal();
  //P_tilde2_sigma ~ std_normal();
  //P_tilde2_raw ~ std_normal();

  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > 1) {
      (ymat[country[i], time[i]] - ymat[country[i], time[i] - 1]) ~ normal(transition_function[country[i], time[i]], epsilon_scale);
    }
  }
}
generated quantities {
  matrix[C, T] eta;
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
    eta[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];

    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition = rate_spline(eta[c, t - 1], P_tilde, P_tilde2, a[c,], ext_knots, num_basis, spline_degree);
      eta[c, t] = eta[c, t - 1] + transition + error;
    }

    for(i in 1:num_grid) {
      transition_function_pred[c, i] = rate_spline(grid[i], 0, 1, a[c,], ext_knots, num_basis, spline_degree);
    }
  }
}
