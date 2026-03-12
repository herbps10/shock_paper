functions {
  vector rate_double_logistic(vector x, vector Delta1, vector Delta2, vector Delta3, vector Delta4, vector k, vector z) {
    real A1 = 4.4;
    real A2 = 0.5;
    return k .* inv(1 + exp(-A1 .* inv(Delta2) .* (x - Delta1 - A2 * Delta2))) + (z - k) .* inv(1 + exp(-A1 * inv(Delta4) .* (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
  }
  
  real shock_rng(real nu_local, real c_slab, real tau) {
    real shock_raw_pred = normal_rng(0, 1);
    real local_shrinkage_pred = student_t_rng(nu_local, 0, 1);
    real truncated_local_shrinkage_pred = sqrt(square(c_slab) * square(local_shrinkage_pred) ./ (square(c_slab) + square(tau) * square(local_shrinkage_pred)));
    return shock_raw_pred * truncated_local_shrinkage_pred * tau;
  } 
  
  real inv_logit_adjustment(real x) {
    return x - 2 * log1p_exp(x);
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

data {
  int C; // Number of countries
  int T; // Number of time points
  int Tpred; // Total number of timepoints
  
  matrix[C, T] y;
  
  int num_grid;
  vector[num_grid] grid;
  
  int<lower=0, upper=1> hierarchical;
  
  real<lower=0> outlier_threshold;
  
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
}
transformed data {
  int n_shocks = C * (T - 1);
  real<lower=1> nu_global = 1;
  real<lower=1> nu_local = 1;
  
  vector[C * (T - 1)] diff = to_vector(y[, 2:T] - y[, 1:(T - 1)]);
  
  int n_below_threshold = 0;
  for(i in 1:(C * (T - 1))) {
    n_below_threshold += (abs(diff[i]) < outlier_threshold) ? 1 : 0;
  }
  array[n_below_threshold] int indices_below_threshold;
  
  {
    int index = 1;
    for(i in 1:(C * (T - 1))) {
      if(abs(diff[i]) < outlier_threshold) {
        indices_below_threshold[index] = i;
        index += 1;
      }
    }
  }
}
parameters {
  //real log_epsilon_scale;
  real<lower=0.01> epsilon_scale;
  
  array[hierarchical] vector[C] raw_Delta1;
  array[hierarchical] vector[C] raw_Delta2;
  array[hierarchical] vector[C] raw_Delta3;
  array[hierarchical] vector[C] raw_Delta4;
  array[hierarchical] vector[C] raw_k;
  array[hierarchical] vector[C] raw_z;
  
  array[hierarchical] real mu_Delta1;
  array[hierarchical] real mu_Delta2;
  array[hierarchical] real mu_Delta3;
  array[hierarchical] real mu_Delta4;
  array[hierarchical] real mu_k;
  array[hierarchical] real mu_z;
  
  array[hierarchical] real<lower=0.01, upper=5> sigma_Delta1;
  array[hierarchical] real<lower=0.01, upper=5> sigma_Delta2;
  array[hierarchical] real<lower=0.01, upper=5> sigma_Delta3;
  array[hierarchical] real<lower=0.01, upper=5> sigma_Delta4;
  array[hierarchical] real<lower=0.01, upper=1> sigma_k;
  array[hierarchical] real<lower=0.01, upper=1> sigma_z;
  
  array[1 - hierarchical] vector<lower=0, upper=100>[C]  constrained_Delta1;
  array[1 - hierarchical] vector<lower=0, upper=100>[C]  constrained_Delta2;
  array[1 - hierarchical] vector<lower=0, upper=100>[C]  constrained_Delta3;
  array[1 - hierarchical] vector<lower=0, upper=100>[C]  constrained_Delta4;
  array[1 - hierarchical] vector<lower=0, upper=10>[C]   constrained_k;
  array[1 - hierarchical] vector<lower=0, upper=1.15>[C] constrained_z;
  
  vector[n_shocks] shock_raw;
  vector<lower=0>[n_shocks] lambda;
  real<lower=0> tau;
  real<lower=0> caux;
}

transformed parameters {
  //real epsilon_scale = exp(log_epsilon_scale);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  
 vector<lower=0, upper=100>[C] Delta1; 
 vector<lower=0, upper=100>[C] Delta2;
 vector<lower=0, upper=100>[C] Delta3;
 vector<lower=0, upper=100>[C] Delta4;
 vector<lower=0, upper=10>[C] k;
 vector<lower=0, upper=1.15>[C] z; 
 
 if(hierarchical) {
    Delta1 = inv_logit(mu_Delta1[1] + sigma_Delta1[1] * raw_Delta1[1]) * 100;
    Delta2 = inv_logit(mu_Delta2[1] + sigma_Delta2[1] * raw_Delta2[1]) * 100;
    Delta3 = inv_logit(mu_Delta3[1] + sigma_Delta3[1] * raw_Delta3[1]) * 100;
    Delta4 = inv_logit(mu_Delta4[1] + sigma_Delta4[1] * raw_Delta4[1]) * 90 + 10;
    k      = inv_logit(mu_k[1] + sigma_k[1] * raw_k[1]) * 10;
    z      = inv_logit(mu_z[1] + sigma_z[1] * raw_z[1]) * 1.15;
  }
  else {
    Delta1 = constrained_Delta1[1];
    Delta2 = constrained_Delta2[1];
    Delta3 = constrained_Delta3[1];
    Delta4 = constrained_Delta4[1];
    k      = constrained_k[1];
    z      = constrained_z[1];
  }
  
  matrix[C, T - 1] shock = rep_matrix(0, C, T - 1);
  real<lower=0> c_slab = slab_scale * sqrt(caux);
  vector<lower=0>[n_shocks] lambda_tilde = sqrt(c_slab^2 * square (lambda) ./ (c_slab^2 + tau^2 * square(lambda)));
  shock = to_matrix(shock_raw .* lambda_tilde * tau, C, T - 1);
  
  transition_function = to_matrix(
    rate_double_logistic(
      to_vector(y[, 1:(T - 1)]),
      rep_vector_times(Delta1, T - 1),
      rep_vector_times(Delta2, T - 1),
      rep_vector_times(Delta3, T - 1),
      rep_vector_times(Delta4, T - 1),
      rep_vector_times(k, T - 1),
      rep_vector_times(z, T - 1)
    ), C, T - 1) + shock;
}

model {
  // here inv gamma is on SD, should be on variance instead
  epsilon_scale ~ inv_gamma(1, 1);
  //epsilon_scale ~ normal(0, 5);
  
  if(hierarchical == 0) {
    Delta1 ~ normal(15.77, 10);
    Delta2 ~ normal(40.97, 10);
    Delta3 ~ normal(0.21, 10);
    Delta4 ~ normal(19.82, 10);
    k ~ normal(2.93, 5);
    z ~ normal(0.4, 0.5);
  }
  else {
    raw_Delta1[1] ~ std_normal();
    raw_Delta2[1] ~ std_normal();
    raw_Delta3[1] ~ std_normal();
    raw_Delta4[1] ~ std_normal();
    raw_k[1]      ~ std_normal();
    raw_z[1]      ~ std_normal();
    
    inv_logit(mu_Delta1[1]) * 100 ~ normal(15.77, 10);
    inv_logit(mu_Delta2[1]) * 100 ~ normal(40.97, 10);
    inv_logit(mu_Delta3[1]) * 100 ~ normal(0.21, 10);
    inv_logit(mu_Delta4[1]) * 90 + 10 ~ normal(19.82, 10);
    inv_logit(mu_k[1]) * 10 ~ normal(2.93, 5);
    inv_logit(mu_z[1]) * 1.15 ~ normal(0.4, 0.5);
    
    target += inv_logit_adjustment(mu_Delta1[1]);
    target += inv_logit_adjustment(mu_Delta2[1]);
    target += inv_logit_adjustment(mu_Delta3[1]);
    target += inv_logit_adjustment(mu_Delta4[1]);
    target += inv_logit_adjustment(mu_k[1]);
    target += inv_logit_adjustment(mu_z[1]);
  }
  
  shock_raw ~ std_normal();
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global);
  
  if(outlier_threshold < 1000) {
    diff[indices_below_threshold] ~ normal(to_vector(transition_function)[indices_below_threshold], sqrt(epsilon_scale));
  }
  else {
    diff ~ normal(to_vector(transition_function), sqrt(epsilon_scale));
  }
}
generated quantities {
  matrix[C, Tpred] eta;
  matrix[C, Tpred] eta_crisisfree;
  matrix[C, Tpred - 1] shock2 = rep_matrix(0, C, Tpred - 1);
  
  matrix[C, num_grid] transition_function_pred;
  
  eta[1:C, 1:T] = y;
  eta_crisisfree[1:C, 1:T] = y;
  
  shock2[1:C, 1:(T - 1)] = shock;
  
  for(t in T:(Tpred)) {
    for(c in 1:C) {
      shock2[c, t - 1] = shock_rng(nu_local, c_slab, tau);
    }
    
    vector[C] transition = rate_double_logistic(eta[, t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
    for(c in 1:C) {
      real error = normal_rng(0, epsilon_scale);
      eta[c, t] = eta[c, t - 1] + transition[c] + error + shock2[c, t - 1];
    }
    
    vector[C] transition_crisisfree = rate_double_logistic(eta_crisisfree[, t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
    for(c in 1:C) {
      real error = normal_rng(0, epsilon_scale);
      eta_crisisfree[c, t] = eta_crisisfree[c, t - 1] + transition_crisisfree[c] + error;
    }
  }
  
  for(i in 1:num_grid) {
    transition_function_pred[, i] = rate_double_logistic(rep_vector(grid[i], C), Delta1, Delta2, Delta3, Delta4, k, z);
  }
}
