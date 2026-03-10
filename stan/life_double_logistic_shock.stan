functions {
  real rate_double_logistic(real x, real Delta1, real Delta2, real Delta3, real Delta4, real k, real z) {
    real A1 = 4.4;
    real A2 = 0.5;
    return k / (1 + exp(-A1 / Delta2 * (x - Delta1 - A2 * Delta2))) + (z - k) / (1 + exp(-A1 / Delta4 * (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
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
  
  vector[N] y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=1> held_out;
  
  int num_grid;
  vector[num_grid] grid;
  
  real<lower=0> outlier_threshold;
  
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
}
transformed data {
  matrix[C, t_last] ymat = rep_matrix(0, C, t_last);
  
  array[C] int final_observed = rep_array(0, C);
  int n_shocks = 0;
  
  real<lower=1> nu_global = 1;
  real<lower=1> nu_local = 1;
  
  for(i in 1:N) {
    ymat[country[i], time[i]] = y[i];
    
    if(held_out[i] == 0 && time[i] > final_observed[country[i]]) {
      final_observed[country[i]] = time[i];
    }
  }
  
  for(c in 1:C) {
    n_shocks += final_observed[c] - 1;
  }
}

parameters {
  real<lower=0> epsilon_scale;
  
  vector<lower=0, upper=100>[C] Delta1;
  vector<lower=0, upper=100>[C] Delta2;
  vector<lower=0, upper=100>[C] Delta3;
  vector<lower=10, upper=100>[C] Delta4;
  vector<lower=0, upper=10>[C] k;
  vector<lower=0, upper=1.15>[C] z;
  
  real <lower=0, upper=100> mu_Delta1;
  real <lower=0, upper=100> mu_Delta2;
  real <lower=0, upper=100> mu_Delta3;
  real <lower=10, upper=100> mu_Delta4;
  real <lower=0, upper=10> mu_k;
  real <lower=0, upper=1.15> mu_z;
  
  real<lower=0, upper=100> sigma_Delta1;
  real<lower=0, upper=100> sigma_Delta2;
  real<lower=0, upper=100> sigma_Delta3;
  real<lower=0, upper=100> sigma_Delta4;
  
  real<lower=0, upper=100> sigma_k;
  real<lower=0, upper=10> sigma_z;
  
  vector<upper=0>[n_shocks] shock_raw;
  real<lower=0> global_shrinkage;
  vector<lower=0>[n_shocks] local_shrinkage; // called lambda in paper
  real<lower=0> caux;
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, t_last] gamma = rep_matrix(0, C, t_last);
  
  matrix[C, t_last] shock = rep_matrix(0, C, t_last);
  
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
  
  for(c in 1:C) {
    for(t in 2:final_observed[c]) {
      transition_function[c, t] = rate_double_logistic(ymat[c, t - 1] - shock[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      gamma[c, t] = transition_function[c, t] + shock[c, t] - shock[c, t - 1];
    }
  }
}

model {
  // here inv gamma is on SD, should be on variance instead
  // epsilon_scale ~ inv_gamma(0.1, 0.1);
  epsilon_scale ~ normal(0, 5);
  
  sigma_Delta1 ~ inv_gamma(0.5, 0.5);
  sigma_Delta2 ~ inv_gamma(0.5, 0.5);
  sigma_Delta3 ~ inv_gamma(0.5, 0.5);
  sigma_Delta4 ~ inv_gamma(0.5, 0.5);
  sigma_k ~ inv_gamma(0.5, 0.5);
  sigma_z ~ inv_gamma(0.5, 0.5);
  
  Delta1 ~ normal(mu_Delta1, sqrt(sigma_Delta1)) T[0, 100];
  Delta2 ~ normal(mu_Delta2, sqrt(sigma_Delta2)) T[0, 100];
  Delta3 ~ normal(mu_Delta3, sqrt(sigma_Delta3)) T[0, 100];
  Delta4 ~ normal(mu_Delta4, sqrt(sigma_Delta4)) T[10, 100];
  k   ~ normal(mu_k, sqrt(sigma_k)) T[0, 10];
  z   ~ normal(mu_z, sqrt(sigma_z)) T[0, 1.15];
  
  mu_Delta1 ~ normal(15.77, 10) T[0, 100];
  mu_Delta2 ~ normal(40.97, 10) T[0, 100];
  mu_Delta3 ~ normal(0.21, 10) T[0, 100];
  mu_Delta4 ~ normal(19.82, 10) T[10, 100];
  mu_k ~ normal(2.93, 5) T[0, 10];
  mu_z ~ normal(0.4, 0.5) T[0, 1.15];
  
  shock_raw ~ std_normal();
  local_shrinkage ~ student_t(nu_local, 0, 1);
  global_shrinkage ~ student_t(nu_global, 0, scale_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > 1) {
      real obs = ymat[country[i], time[i]] - ymat[country[i], time[i] - 1];
      if(abs(obs) < outlier_threshold) {
        obs ~ normal(gamma[country[i], time[i]], epsilon_scale);
      }
    }
  }
}
generated quantities {
  matrix[C, T] eta;
  matrix[C, T] eta_crisisfree;
  matrix[C, T] shock2;
  matrix[C, num_grid] transition_function_pred;
  
  vector[num_grid] transition_function_mean = rep_vector(0, num_grid);
  //for(i in 1:num_grid) {
  //  transition_function_mean[i] = rate_spline(grid[i], 0, 1, to_row_vector(a_mean), ext_knots, num_basis, spline_degree);
  //}
  
  for(c in 1:C) {
    eta_crisisfree[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]] - shock[c, 1:final_observed[c]];
    
    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition_crisisfree = rate_double_logistic(eta_crisisfree[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      eta_crisisfree[c, t] = eta_crisisfree[c, t - 1] + transition_crisisfree + error;
    }
    
    eta[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];
    shock2[c, 1:(final_observed[c])] = shock[c, 1:(final_observed[c])];
    
    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      shock2[c, t] = shock_rng(nu_local, c_slab, global_shrinkage);
      real transition = rate_double_logistic(eta[c, t - 1] - shock2[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      eta[c, t] = eta[c, t - 1] + transition + error + shock2[c, t] - shock2[c, t - 1];
    }
    
    for(i in 1:num_grid) {
      transition_function_pred[c, i] = rate_double_logistic(grid[i], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
    }
  }
}
