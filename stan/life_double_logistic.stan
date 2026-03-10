functions {
  real rate_double_logistic(real x, real Delta1, real Delta2, real Delta3, real Delta4, real k, real z) {
    real A1 = 4.4;
    real A2 = 0.5;
    return k / (1 + exp(-A1 / Delta2 * (x - Delta1 - A2 * Delta2))) + (z - k) / (1 + exp(-A1 / Delta4 * (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
  }
  
  real inv_logit_adjustment(real x) {
    return x - 2 * log(exp(x) + 1);
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int C; // Number of countries
  int t_last;
  
  vector[N] y;                             // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=1> held_out;
  
  int num_grid;
  vector[num_grid] grid;
  
  real<lower=0> outlier_threshold;
}
transformed data {
  matrix[C, t_last] ymat = rep_matrix(0, C, t_last);
  
  array[C] int final_observed = rep_array(0, C);
  
  for(i in 1:N) {
    ymat[country[i], time[i]] = y[i];
    
    if(held_out[i] == 0 && time[i] > final_observed[country[i]]) {
      final_observed[country[i]] = time[i];
    }
  }
}

parameters {
  real<lower=0.5,upper=1> epsilon_scale;
  
  vector<lower=-5, upper=5>[C] raw_Delta1;
  vector<lower=-5, upper=5>[C] raw_Delta2;
  vector<lower=-5, upper=5>[C] raw_Delta3;
  vector<lower=-5, upper=5>[C] raw_Delta4;
  vector<lower=-5, upper=5>[C] raw_k;
  vector<lower=-5, upper=5>[C] raw_z;
  
  real<lower=-5, upper=5> mu_Delta1;
  real<lower=-5, upper=5> mu_Delta2;
  real<lower=-5, upper=5> mu_Delta3;
  real<lower=-5, upper=5> mu_Delta4;
  real<lower=-5, upper=5> mu_k;
  real<lower=-5, upper=5> mu_z;
  
  real<lower=0.1, upper=5> sigma_Delta1;
  real<lower=0.1, upper=5> sigma_Delta2;
  real<lower=0.1, upper=5> sigma_Delta3;
  real<lower=0.1, upper=5> sigma_Delta4;
  real<lower=0.1, upper=1> sigma_k;
  real<lower=0.1, upper=1> sigma_z;
}

transformed parameters {
  matrix[C, t_last] transition_function = rep_matrix(0, C, t_last);
  matrix[C, t_last] gamma = rep_matrix(0, C, t_last);
  
  vector<lower=0, upper=100>[C]  Delta1 = inv_logit(mu_Delta1 + sigma_Delta1 * raw_Delta1) * 100;
  vector<lower=0, upper=100>[C]  Delta2 = inv_logit(mu_Delta2 + sigma_Delta2 * raw_Delta2) * 100;
  vector<lower=0, upper=100>[C]  Delta3 = inv_logit(mu_Delta3 + sigma_Delta3 * raw_Delta3) * 100;
  vector<lower=10, upper=100>[C] Delta4 = inv_logit(mu_Delta4 + sigma_Delta4 * raw_Delta4) * 90 + 10;
  vector<lower=0, upper=10>[C]   k      = inv_logit(mu_k + sigma_k * raw_k) * 10;
  vector<lower=0, upper=1.15>[C] z      = inv_logit(mu_z + sigma_z * raw_z) * 1.15;
  
  for(c in 1:C) {
    for(t in 2:final_observed[c]) {
      transition_function[c, t] = rate_double_logistic(ymat[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      gamma[c, t] = transition_function[c, t];
    }
  }
}

model {
  // here inv gamma is on SD, should be on variance instead
  // epsilon_scale ~ inv_gamma(0.1, 0.1);
  epsilon_scale ~ normal(0, 5);
  
  raw_Delta1 ~ std_normal();
  raw_Delta2 ~ std_normal();
  raw_Delta3 ~ std_normal();
  raw_Delta4 ~ std_normal();
  raw_k      ~ std_normal();
  raw_z      ~ std_normal();
  
  inv_logit(mu_Delta1) * 100 ~ normal(15.77, 10) T[0, 100];
  inv_logit(mu_Delta2) * 100 ~ normal(40.97, 10) T[0, 100];
  inv_logit(mu_Delta3) * 100 ~ normal(0.21, 10) T[0, 100];
  inv_logit(mu_Delta4) * 90 + 10 ~ normal(19.82, 10) T[10, 100];
  inv_logit(mu_k) * 10 ~ normal(2.93, 5) T[0, 10];
  inv_logit(mu_z) * 1.15 ~ normal(0.4, 0.5) T[0, 1.15];
  
  target += inv_logit_adjustment(mu_Delta1);
  target += inv_logit_adjustment(mu_Delta2);
  target += inv_logit_adjustment(mu_Delta3);
  target += inv_logit_adjustment(mu_Delta4);
  target += inv_logit_adjustment(mu_k);
  target += inv_logit_adjustment(mu_z);
  
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
  matrix[C, num_grid] transition_function_pred;
  vector[num_grid] transition_function_mean = rep_vector(0, num_grid);

  //for(i in 1:num_grid) {
  //  transition_function_mean[i] = rate_spline(grid[i], 0, 1, to_row_vector(a_mean), ext_knots, num_basis, spline_degree);
  //}
  
  for(c in 1:C) {
    eta_crisisfree[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];
    
    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition_crisisfree = rate_double_logistic(eta_crisisfree[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      eta_crisisfree[c, t] = eta_crisisfree[c, t - 1] + transition_crisisfree + error;
    }
    
    eta[c, 1:final_observed[c]] = ymat[c, 1:final_observed[c]];
    
    for(t in (final_observed[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition = rate_double_logistic(eta[c, t - 1], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
      eta[c, t] = eta[c, t - 1] + transition + error;
    }
    
    for(i in 1:num_grid) {
      transition_function_pred[c, i] = rate_double_logistic(grid[i], Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c]);
    }
  }
}
