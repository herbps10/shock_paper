functions {
  vector rate_double_logistic(vector x, vector Delta1, vector Delta2, vector Delta3, vector Delta4, vector k, vector z) {
    real A1 = 4.4;
    real A2 = 0.5;
    return k .* inv(1 + exp(-A1 .* inv(Delta2) .* (x - Delta1 - A2 * Delta2))) + (z - k) .* inv(1 + exp(-A1 * inv(Delta4) .* (x - Delta1 - Delta2 - Delta3 - A2 * Delta4)));
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
  int<lower=0, upper=1> centered;
  
  real<lower=0> outlier_threshold;
}
transformed data {
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

  real Delta1_lower = 0;  real Delta1_upper = 25; real Delta1_range = Delta1_upper - Delta1_lower;
  real Delta2_lower = 25; real Delta2_upper = 50; real Delta2_range = Delta2_upper - Delta2_lower;
  real Delta3_lower = 0;  real Delta3_upper = 10; real Delta3_range = Delta3_upper - Delta3_lower;
  real Delta4_lower = 5;  real Delta4_upper = 30; real Delta4_range = Delta4_upper - Delta4_lower;
  real k_lower = 0; real k_upper = 10; real k_range = k_upper - k_lower;
  real z_lower = 0; real z_upper = 1.15/5.0; real z_range = z_upper - z_lower;
  
  real prior_mu_Delta1 = logit((15.77 - Delta1_lower) / Delta1_range);
  real prior_mu_Delta2 = logit((40.97 - Delta2_lower) / Delta1_range);
  real prior_mu_Delta3 = logit(( 0.21 - Delta3_lower) / Delta3_range);
  real prior_mu_Delta4 = logit((19.82 - Delta4_lower) / Delta4_range);
  real prior_mu_k      = logit((2.93 - k_lower) / k_range);
  real prior_mu_z      = logit(( 0.4/5.0 - z_lower) / z_range);
}
parameters {
  real<lower=0> epsilon_sigma;
  
  array[hierarchical * (1 - centered)] vector[C] raw_Delta1;
  array[hierarchical * (1 - centered)] vector[C] raw_Delta2;
  array[hierarchical * (1 - centered)] vector[C] raw_Delta3;
  array[hierarchical * (1 - centered)] vector[C] raw_Delta4;
  array[hierarchical * (1 - centered)] vector[C] raw_k;
  array[hierarchical * (1 - centered)] vector[C] raw_z;
  
  array[hierarchical * centered] vector[C] Delta1_logit;
  array[hierarchical * centered] vector[C] Delta2_logit;
  array[hierarchical * centered] vector[C] Delta3_logit;
  array[hierarchical * centered] vector[C] Delta4_logit;
  array[hierarchical * centered] vector[C] k_logit;
  array[hierarchical * centered] vector[C] z_logit;
  
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

  array[1 - hierarchical] vector<lower=Delta1_lower, upper=Delta1_upper>[C] constrained_Delta1;
  array[1 - hierarchical] vector<lower=Delta2_lower, upper=Delta2_upper>[C] constrained_Delta2;
  array[1 - hierarchical] vector<lower=Delta3_lower, upper=Delta3_upper>[C] constrained_Delta3;
  array[1 - hierarchical] vector<lower=Delta4_lower, upper=Delta4_upper>[C] constrained_Delta4;
  array[1 - hierarchical] vector<lower=k_lower, upper=k_upper>[C] constrained_k;
  array[1 - hierarchical] vector<lower=z_lower, upper=z_upper>[C] constrained_z;
}

transformed parameters {
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  
  vector[C] Delta1; 
  vector[C] Delta2;
  vector[C] Delta3;
  vector[C] Delta4;
  vector[C] k;
  vector[C] z; 
  
  if(hierarchical) {
    if(centered) {
      Delta1 = inv_logit(Delta1_logit[1]) * Delta1_range + Delta1_lower;
      Delta2 = inv_logit(Delta2_logit[1]) * Delta2_range + Delta2_lower;
      Delta3 = inv_logit(Delta3_logit[1]) * Delta3_range + Delta3_lower;
      Delta4 = inv_logit(Delta4_logit[1]) * Delta4_range + Delta4_lower;
      k      = inv_logit(k_logit[1]) * k_range + k_lower;
      z      = inv_logit(z_logit[1]) * z_range + z_lower;
    }
    else {
      Delta1 = inv_logit(mu_Delta1[1] + sigma_Delta1[1] * raw_Delta1[1]) * Delta1_range + Delta1_lower; 
      Delta2 = inv_logit(mu_Delta2[1] + sigma_Delta2[1] * raw_Delta2[1]) * Delta2_range + Delta2_lower; 
      Delta3 = inv_logit(mu_Delta3[1] + sigma_Delta3[1] * raw_Delta3[1]) * Delta3_range + Delta3_lower; 
      Delta4 = inv_logit(mu_Delta4[1] + sigma_Delta4[1] * raw_Delta4[1]) * Delta4_range + Delta4_lower; 
      k      = inv_logit(mu_k[1] + sigma_k[1] * raw_k[1]) * z_range + z_lower;
      z      = inv_logit(mu_z[1] + sigma_z[1] * raw_z[1]) * k_range + k_lower;
    }
  }
  else {
    Delta1 = constrained_Delta1[1];
    Delta2 = constrained_Delta2[1];
    Delta3 = constrained_Delta3[1];
    Delta4 = constrained_Delta4[1];
    k      = constrained_k[1];
    z      = constrained_z[1];
  }
  
  transition_function = to_matrix(
    rate_double_logistic(
      to_vector(y[, 1:(T - 1)]),
      rep_vector_times(Delta1, T - 1),
      rep_vector_times(Delta2, T - 1),
      rep_vector_times(Delta3, T - 1),
      rep_vector_times(Delta4, T - 1),
      rep_vector_times(k, T - 1),
      rep_vector_times(z, T - 1)
    ), C, T - 1);
}

model {
  // here inv gamma is on SD, should be on variance instead
  epsilon_sigma ~ normal(0, 2);
  //epsilon_sigma ~ normal(0, 5);
  
  if(hierarchical == 0) {
    Delta1 ~ normal(15.77, 10);
    Delta2 ~ normal(40.97, 10);
    Delta3 ~ normal(0.21, 10);
    Delta4 ~ normal(19.82, 10);
    k ~ normal(2.93, 5);
    z ~ normal(0.4/5.0, 0.5);
  }
  else {
    mu_Delta1[1] ~ normal(prior_mu_Delta1, 2);
    mu_Delta2[1] ~ normal(prior_mu_Delta2, 2);
    mu_Delta3[1] ~ normal(prior_mu_Delta3, 2);
    mu_Delta4[1] ~ normal(prior_mu_Delta4, 2);
    mu_k[1] ~ normal(prior_mu_k, 2);
    mu_z[1] ~ normal(prior_mu_z, 2);
    
    sigma_Delta1[1] ~ normal(0, 2);
    sigma_Delta2[1] ~ normal(0, 2);
    sigma_Delta3[1] ~ normal(0, 2);
    sigma_Delta4[1] ~ normal(0, 2);
    sigma_k[1] ~ normal(0, 1);
    sigma_z[1] ~ normal(0, 1);
    
    if(centered == 1) {
      Delta1_logit[1] ~ normal(mu_Delta1[1], sigma_Delta1[1]);
      Delta2_logit[1] ~ normal(mu_Delta2[1], sigma_Delta2[1]);
      Delta3_logit[1] ~ normal(mu_Delta3[1], sigma_Delta3[1]);
      Delta4_logit[1] ~ normal(mu_Delta4[1], sigma_Delta4[1]);
      k_logit[1] ~ normal(k_logit[1], sigma_k[1]);
      z_logit[1] ~ normal(z_logit[1], sigma_z[1]);
    }
    else {
      raw_Delta1[1] ~ std_normal();
      raw_Delta2[1] ~ std_normal();
      raw_Delta3[1] ~ std_normal();
      raw_Delta4[1] ~ std_normal();
      raw_k[1]      ~ std_normal();
      raw_z[1]      ~ std_normal();
    }
  }
  
  if(outlier_threshold < 1000) {
    diff[indices_below_threshold] ~ normal(to_vector(transition_function)[indices_below_threshold], epsilon_sigma);
  }
  else {
    diff ~ normal(to_vector(transition_function), epsilon_sigma);
  }
}
generated quantities {
  matrix[C, Tpred] eta;
  matrix[C, num_grid] transition_function_pred;
  
  eta[1:C, 1:T] = y;
  
  for(t in T:Tpred) {
    vector[C] transition = rate_double_logistic(eta[, t - 1], Delta1, Delta2, Delta3, Delta4, k, z);
    for(c in 1:C) {
      real error = normal_rng(0, sqrt(epsilon_variance));
      eta[c, t] = eta[c, t - 1] + transition[c] + error;
    }
  }
  
  for(i in 1:num_grid) {
    transition_function_pred[, i] = rate_double_logistic(rep_vector(grid[i], C), Delta1, Delta2, Delta3, Delta4, k, z);
  }
}
