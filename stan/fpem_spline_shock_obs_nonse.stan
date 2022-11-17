functions {
#include ./fill_AR.stan
#include ./scale_blocks.stan
#include ./deboor.stan

  // From the Stan User's Guide 2.28: https://mc-stan.org/docs/2_28/stan-users-guide/truncated-random-number-generation.html.
  // Released under CC-BY 4.0 license (https://creativecommons.org/licenses/by/4.0/legalcode)
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }

  // From the Stan User's Guide 2.28: https://mc-stan.org/docs/2_28/stan-users-guide/truncated-random-number-generation.html.
  // Released under CC-BY 4.0 license (https://creativecommons.org/licenses/by/4.0/legalcode)
  real normal_lub_cdf(real x, real mu, real sigma, real lb, real ub) {
    real p_x  = normal_cdf(x  | mu, sigma);
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    return (p_x - p_lb) / p_ub;
  }

  real rate_spline(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
    return deboor(P / P_tilde, ext_knots, a, spline_degree);
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int C; // Number of countries
  int S; // Number of sources

  array[C] int<lower=1, upper=T> t_star;
  int<lower=1, upper=T> t_last;

  array[N] real y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=S> source;   // Source of each observation
  array[N] int<lower=0, upper=1> held_out; // Whether to hold out each observation

  array[N] real<lower=0> s;                // Standard deviation

  int num_knots;
  vector[num_knots] knots;

  int spline_degree;

  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;

  // hierarchical P_tilde
  int<lower=0> P_tilde_n_terms;
  int<lower=0> P_tilde_n_re;
  array[P_tilde_n_re] int<lower=1, upper=P_tilde_n_terms> P_tilde_re_start;
  array[P_tilde_n_re] int<lower=1, upper=P_tilde_n_terms> P_tilde_re_end;
  matrix[C, P_tilde_n_terms] P_tilde_model_matrix;

  // hierarchical Omega;
  int<lower=0> Omega_n_terms;
  int<lower=0> Omega_n_re;
  array[Omega_n_re] int<lower=1, upper=Omega_n_terms> Omega_re_start;
  array[Omega_n_re] int<lower=1, upper=Omega_n_terms> Omega_re_end;
  matrix[C, Omega_n_terms] Omega_model_matrix;


  // hierarchical a
  int<lower=0> a_n_terms;
  int<lower=0> a_n_re;
  array[a_n_re] int<lower=1, upper=a_n_terms> a_re_start;
  array[a_n_re] int<lower=1, upper=a_n_terms> a_re_end;
  matrix[C, a_n_terms] a_model_matrix;

  real a_lower_bound;
  real a_upper_bound;

  int<lower=0, upper=1> smoothing; // 0: no smoothing component. 1: AR(1) smoothing component.
  
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
}
transformed data {
  vector[rows(csr_extract_w(P_tilde_model_matrix))] P_tilde_model_matrix_w    = csr_extract_w(P_tilde_model_matrix);
  array[size(csr_extract_v(P_tilde_model_matrix))] int P_tilde_model_matrix_v = csr_extract_v(P_tilde_model_matrix);
  array[size(csr_extract_u(P_tilde_model_matrix))] int P_tilde_model_matrix_u = csr_extract_u(P_tilde_model_matrix);

  vector[rows(csr_extract_w(Omega_model_matrix))] Omega_model_matrix_w     = csr_extract_w(Omega_model_matrix);
  array[size(csr_extract_v(Omega_model_matrix))] int Omega_model_matrix_v  = csr_extract_v(Omega_model_matrix);
  array[size(csr_extract_u(Omega_model_matrix))] int Omega_model_matrix_u  = csr_extract_u(Omega_model_matrix);

  vector[rows(csr_extract_w(a_model_matrix))] a_model_matrix_w             = csr_extract_w(a_model_matrix);
  array[size(csr_extract_v(a_model_matrix))] int a_model_matrix_v          = csr_extract_v(a_model_matrix);
  array[size(csr_extract_u(a_model_matrix))] int a_model_matrix_u          = csr_extract_u(a_model_matrix);


  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  int num_constrained_zero = spline_degree + 1;
  
  real<lower=1> nu_global = 1;
  real<lower=1> nu_local = 1;
  real<lower=0> slab_df = 1;

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;
}

parameters {
  // Spline coefficients
  matrix[a_n_terms, num_basis - num_constrained_zero] a_raw;
  matrix<lower=0>[a_n_re - 1, num_basis - num_constrained_zero] a_sigma_raw;

  // P_tilde
  vector[P_tilde_n_terms] P_tilde_raw;
  vector<lower=0>[P_tilde_n_re - 1] P_tilde_sigma_raw;

  // Omega
  vector[Omega_n_terms] Omega_raw;
  vector<lower=0>[Omega_n_re - 1] Omega_sigma_raw;

  // Smoothing component
  matrix[C * smoothing, T * smoothing] epsilon_innovation;
  array[smoothing] real<lower=0, upper=1> est_rho;
  array[smoothing] real<lower=0> est_tau;

  // Data model
  array[S] real<lower=0> nonse;
  
  vector<upper=0>[C * t_last] shock_raw;
  real<lower=0> global_shrinkage;
  vector<lower=0>[C * t_last] local_shrinkage; // called lambda in paper
  real<lower=0> caux;
  
  vector<lower=0>[N] obs_nonse_raw;
  real<lower=0> obs_global_shrinkage;
  vector<lower=0>[N] obs_local_shrinkage; // called lambda in paper
  real<lower=0> obs_caux;
}

transformed parameters {
  matrix[C, T] eta;
  matrix[C, num_basis] a;
  array[smoothing] real rho;
  array[smoothing] real tau;
  matrix[C * smoothing, T * smoothing] epsilon;

  vector<lower=0, upper=1>[C] P_tilde;
  vector[P_tilde_n_terms] P_tilde_star;

  vector[C] Omega;
  vector[Omega_n_terms] Omega_star;

  matrix[a_n_terms, num_basis - num_constrained_zero] a_star;

  matrix<lower=0>[a_n_re, num_basis - num_constrained_zero] a_sigma;
  vector<lower=0>[P_tilde_n_re] P_tilde_sigma;
  vector<lower=0>[Omega_n_re] Omega_sigma;

  vector<lower=0>[N] scale;
  
  vector<lower=0>[C * t_last] truncated_local_shrinkage; // called lambda_tilde in paper
  real<lower=0> c_slab;
  matrix[C, t_last] shock;
  
  vector<lower=0>[N] obs_truncated_local_shrinkage; // called lambda_tilde in paper
  real<lower=0> obs_c_slab;
  vector<lower=0>[N] obs_nonse;
  
  
  c_slab = slab_scale * sqrt(caux);
  truncated_local_shrinkage = sqrt(c_slab^2 * square(local_shrinkage) ./ (c_slab^2 + global_shrinkage^2 * square(local_shrinkage)));
  shock = to_matrix(shock_raw .* truncated_local_shrinkage * global_shrinkage, C, t_last);
  
  obs_c_slab = slab_scale * sqrt(obs_caux);
  obs_truncated_local_shrinkage = sqrt(obs_c_slab^2 * square(obs_local_shrinkage) ./ (obs_c_slab^2 + obs_global_shrinkage^2 * square(obs_local_shrinkage)));
  obs_nonse = obs_nonse_raw .* obs_truncated_local_shrinkage * obs_global_shrinkage;

  a_sigma[1, ] = rep_row_vector(1, num_basis - num_constrained_zero);
  if(a_n_re > 1) a_sigma[2:a_n_re, ] = a_sigma_raw;

  Omega_sigma[1] = 1;
  if(Omega_n_re > 1) Omega_sigma[2:Omega_n_re] = Omega_sigma_raw;

  P_tilde_sigma[1] = 1;
  if(P_tilde_n_re > 1) P_tilde_sigma[2:P_tilde_n_re] = P_tilde_sigma_raw;

  if(smoothing == 1) {
    rho[1] = est_rho[1];
    tau[1] = est_tau[1];
  }

  // P_tilde
  P_tilde_star = scale_blocks(P_tilde_raw, P_tilde_sigma, P_tilde_re_start, P_tilde_re_end);
  P_tilde = 0.5 + 0.45 * inv_logit(csr_matrix_times_vector(C, P_tilde_n_terms, P_tilde_model_matrix_w, P_tilde_model_matrix_v, P_tilde_model_matrix_u, P_tilde_star));

  // Omega
  Omega_star = scale_blocks(Omega_raw, Omega_sigma, Omega_re_start, Omega_re_end);
  Omega = csr_matrix_times_vector(C, Omega_n_terms, Omega_model_matrix_w, Omega_model_matrix_v, Omega_model_matrix_u, Omega_star);

  // Initialize the non-zero spline coefficients
  for(i in 1:(num_basis - num_constrained_zero)) {
    a_star[, i] = scale_blocks(a_raw[, i], a_sigma[, i], a_re_start, a_re_end);
    a[, i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(csr_matrix_times_vector(C, a_n_terms, a_model_matrix_w, a_model_matrix_v, a_model_matrix_u, a_star[, i]));
  }


  for(c in 1:C) {
    row_vector[T] logit_eta;
    real transition_function;

    a[c, (num_basis - num_constrained_zero + 1):num_basis] = rep_row_vector(0, num_constrained_zero);

    if(smoothing == 1) {
      epsilon[c,] = fill_AR(epsilon_innovation[c, ], rho[1], tau[1], t_star[c]);
    }

    logit_eta[t_star[c]] = Omega[c];

    // Additive smoothing
    if(smoothing == 1) {
      for(t in (t_star[c] + 1):T) {
        transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t - 1] + transition_function + epsilon[c, t];
        if(t <= t_last) logit_eta[t] += shock[c, t];
      }

      for(q in 1:(t_star[c] - 1)) {
        int t = t_star[c] - q;
        transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t + 1] - transition_function - epsilon[c, t + 1];
        if(t <= t_last) logit_eta[t] -= shock[c, t];
      }
    }
    // No smoothing
    else if(smoothing == 0) {
      for(t in (t_star[c] + 1):T) {
        transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t - 1] + transition_function;
        if(t <= t_last) logit_eta[t] += shock[c,t];
      }

      for(q in 1:(t_star[c] - 1)) {
        int t = t_star[c] - q;
        transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t + 1] - transition_function;
        if(t <= t_last) logit_eta[t] -= shock[c,t];
      }
    }

    eta[c, ] = inv_logit(logit_eta);
  }
  for(i in 1:N) {
    scale[i] = sqrt(square(s[i]) + square(nonse[source[i]]) + square(obs_nonse[i]));
  }
}

model {
  // P_tilde
  P_tilde_sigma ~ std_normal();
  P_tilde_raw   ~ std_normal();

  // Omega
  Omega_sigma ~ std_normal();
  Omega_raw   ~ std_normal();


  if(smoothing == 1) {
    est_rho[1] ~ {{RHO_PRIOR}};
    est_tau[1] ~ {{TAU_PRIOR}};
    to_vector(epsilon_innovation) ~ std_normal();
  }

  to_vector(a_raw) ~ std_normal();
  to_vector(a_sigma) ~ std_normal();

  nonse ~ normal(0, 0.1);
  
  shock_raw ~ std_normal();
  local_shrinkage ~ student_t(nu_local, 0, 1);
  global_shrinkage ~ student_t(nu_global, 0, scale_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  obs_nonse_raw ~ std_normal();
  obs_local_shrinkage ~ student_t(nu_local, 0, 1);
  obs_global_shrinkage ~ student_t(nu_global, 0, scale_global);
  obs_caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  for(i in 1:N) {
    if(held_out[i] == 0) {
      y[i] ~ normal(eta[country[i], time[i]], scale[i]) T[0, 1];
    }
  }
}
generated quantities {
  vector[N] y_pred;
  vector[N] resid;
  vector[N] pit;

  for(i in 1:N) {
    y_pred[i]  = normal_lub_rng(eta[country[i], time[i]], scale[i], 0, 1);
    resid[i]   = y_pred[i] - y[i];
    pit[i]     = normal_lub_cdf(y[i] | eta[country[i], time[i]], scale[i], 0, 1);
  }
}
