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

  // hierarchical a
  int<lower=0> a_n_terms;
  int<lower=0> a_n_re;
  int<lower=1, upper=a_n_terms> a_re_start[a_n_re];
  int<lower=1, upper=a_n_terms> a_re_end[a_n_re];
  matrix[C, a_n_terms] a_model_matrix;

  real a_lower_bound;
  real a_upper_bound;

  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;

  int<lower=1> R;
}
transformed data {
  vector[rows(csr_extract_w(a_model_matrix))] a_model_matrix_w             = csr_extract_w(a_model_matrix);
  int a_model_matrix_v[size(csr_extract_v(a_model_matrix))]                = csr_extract_v(a_model_matrix);
  int a_model_matrix_u[size(csr_extract_u(a_model_matrix))]                = csr_extract_u(a_model_matrix);

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
    n_shocks += final_observed[c] - 2;
  }
}

parameters {
  // Spline rate vs. level function
  matrix[a_n_terms, num_basis - num_constrained_zero] a_raw;
  matrix<lower=0>[a_n_re - 1, num_basis - num_constrained_zero] a_sigma_raw;

  //vector[C] Omega_raw;
  real P_tilde2_mu;
  vector[C] P_tilde2_raw;

  real<lower=0> epsilon_scale;

  vector<upper=0>[n_shocks] shock_raw;
  real<lower=0> global_shrinkage_raw;
  //vector<lower=0>[C] global_shrinkage;
  vector<lower=0>[n_shocks] local_shrinkage; // called lambda in paper
  real<lower=0> caux;
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, t_last] gamma = rep_matrix(0, C, t_last);
  matrix[C, num_basis] a;
  //real<lower=0> epsilon_scale = 0.3;

  vector<lower=0>[C] global_shrinkage = rep_vector(global_shrinkage_raw, C);

  //vector[C] Omega = 15 + inv_logit(Omega_raw) * 75;
  vector[C] P_tilde = rep_vector(15, C);
  vector[C] P_tilde2 = inv_logit(P_tilde2_mu + P_tilde2_raw) * 100;

  matrix[a_n_terms, num_basis - num_constrained_zero] a_star;

  vector<lower=0>[n_shocks] truncated_local_shrinkage; // called lambda_tilde in paper
  real<lower=0> c_slab;
  matrix[C, t_last - 1] shock = rep_matrix(0, C, t_last - 1);

  matrix<lower=0>[a_n_re, num_basis - num_constrained_zero] a_sigma;
  a_sigma[1, ] = rep_row_vector(5, num_basis - num_constrained_zero);
  a_sigma[2:a_n_re, ] = a_sigma_raw;

  c_slab = slab_scale * sqrt(caux);
  //truncated_local_shrinkage = sqrt(c_slab^2 * square(local_shrinkage) ./ (c_slab^2 + global_shrinkage^2 * square(local_shrinkage)));
  //shock = to_matrix(shock_raw .* truncated_local_shrinkage * global_shrinkage, C, t_last - 1);

  {
    vector[n_shocks] shock_shrinkage;
    int index = 1;

    for(c in 1:C) {
      truncated_local_shrinkage[index:index + final_observed[c] - 3] = sqrt(c_slab^2 * square(local_shrinkage[index:(index + final_observed[c] - 3)]) ./ (c_slab^2 + global_shrinkage[c]^2 * square(local_shrinkage[index:(index + final_observed[c] - 3)])));
      shock_shrinkage[index:(index + final_observed[c] - 3)] = shock_raw[index:(index + final_observed[c] - 3)] .* truncated_local_shrinkage[index:(index + final_observed[c] - 3)] * global_shrinkage[c];

      shock[c, 1:(final_observed[c] - 2)] = to_row_vector(shock_shrinkage[index:(index + final_observed[c] - 3)]);
      index += (final_observed[c] - 2);
    }
  }

  // Initialize the non-zero spline coefficients
  for(i in 1:(num_basis - num_constrained_zero)) {
    a_star[, i] = scale_blocks(a_raw[, i], a_sigma[, i], a_re_start, a_re_end);
    a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(csr_matrix_times_vector(C, a_n_terms, a_model_matrix_w, a_model_matrix_v, a_model_matrix_u, a_star[, i]));
  }

  for(c in 1:C) {
    a[c, (num_basis - num_constrained_zero + 1):num_basis] = rep_row_vector(0, num_constrained_zero);
    for(t in 2:final_observed[c]) {
      transition_function[c, t] = rate_spline(ymat[c, t - 1], P_tilde[c], P_tilde2[c], a[c,], ext_knots, num_basis, spline_degree);
      gamma[c, t] = transition_function[c, t] + shock[c, t - 1];
    }
  }
}

model {
  to_vector(a_raw) ~ std_normal();
  to_vector(a_sigma_raw) ~ normal(0, 5);

  epsilon_scale ~ normal(0, 1);

  P_tilde2_mu ~ std_normal();
  P_tilde2_raw ~ std_normal();

  shock_raw ~ std_normal();
  local_shrinkage ~ student_t(nu_local, 0, 1);
  global_shrinkage_raw ~ student_t(nu_global, 0, scale_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > 1) {
      (ymat[country[i], time[i]] - ymat[country[i], time[i] - 1]) ~ normal(gamma[country[i], time[i]], epsilon_scale);
    }
  }
}
generated quantities {
  array[C, T, R] real eta;
  array[C, T, 1] real eta_crisisfree;

  for(c in 1:C) {
    eta_crisisfree[c, 1:final_observed[c], 1] = to_array_1d(ymat[c, 1:final_observed[c]]);

    for(t in (final_observed[c] + 1):T) {
      real transition_crisisfree = rate_spline(eta_crisisfree[c, t - 1, 1], P_tilde[c], P_tilde2[c], a[c,], ext_knots, num_basis, spline_degree);
      real error = normal_rng(0, epsilon_scale);

      eta_crisisfree[c, t, 1] = eta_crisisfree[c, t - 1, 1] + transition_crisisfree + error;
    }

    for(r in 1:R) {
      eta[c, 1:final_observed[c], r] = to_array_1d(ymat[c, 1:final_observed[c]]);

      for(t in (final_observed[c] + 1):T) {
        real error = normal_rng(0, epsilon_scale);
        real transition = rate_spline(eta[c, t - 1, r], P_tilde[c], P_tilde2[c], a[c,], ext_knots, num_basis, spline_degree);

        real shock_raw_pred = normal_lub_rng(0, 1, negative_infinity(), 0);
        real local_shrinkage_pred = student_t_rng(nu_local, 0, 1);
        real truncated_local_shrinkage_pred = sqrt(c_slab^2 * square(local_shrinkage_pred) ./ (c_slab^2 + global_shrinkage[c]^2 * square(local_shrinkage_pred)));
        real shock_shrinkage_pred = shock_raw_pred * truncated_local_shrinkage_pred * global_shrinkage[c];

        eta[c, t, r] = eta[c, t - 1, r] + shock_shrinkage_pred + transition + error;
      }
    }

  }
}
