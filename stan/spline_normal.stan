functions {
  real error_rng(vector x, vector epsilon_params) {
    return normal_rng(0, epsilon_params[1]);
  }
  
  real deboor(real x, vector ext_knots, row_vector a, int degree) {
    int k = degree + 1;
    row_vector[degree + 1] d;
    int n_ext_knots = rows(ext_knots);
    
    if (x <= ext_knots[1]) 
      return a[1];
    if (x >= ext_knots[n_ext_knots]) 
      return a[cols(a)];
    
    while (!(ext_knots[k + 1] > x) && k < n_ext_knots - degree - 1) {
      k = k + 1;
    }
    
    d = a[(k - degree) : (k)];
    
    for (r in 2 : (degree + 1)) {
      for (j in (k + r - degree - 1) : k) {
        int j2 = (k + r - degree - 1) - j + k;
        real alpha = (x - ext_knots[j2])
                     / (ext_knots[j2 + 1 + degree - (r - 1)] - ext_knots[j2]);
        
        d[j2 + 1 - k + degree] = (1 - alpha) * d[j2 - k + degree]
                                 + alpha * d[j2 + 1 - k + degree];
      }
    }
    
    return d[degree + 1];
  }
  
  real rate_spline(real P, real P_tilde, real P_tilde2, row_vector alpha,
                   vector ext_knots, int num_basis, int spline_degree) {
    return deboor((P - P_tilde) / (P_tilde2 - P_tilde), ext_knots, alpha,
                  spline_degree);
  }
}
data {
  int C;
  int T;
  int Tpred;
  
  matrix[C, T] y;
  
  int num_grid;
  vector[num_grid] grid;
  
  real<lower=0> outlier_threshold;
  
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  matrix[num_knots + spline_degree - 1, num_grid] B;
}
transformed data {
  vector[C * (T - 1)] diff = to_vector(y[ : , 2 : T] - y[ : , 1 : (T - 1)]);
  
  array[1] int epsilon_params_int;
  epsilon_params_int[1] = 0;
  
  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  
  ext_knots[1 : spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1) : (num_knots + 2 * spline_degree)] = rep_vector(
                                                                    knots[num_knots],
                                                                    spline_degree);
  ext_knots[(spline_degree + 1) : (num_knots + spline_degree)] = knots;
  
  real P_tilde = 5;
  real P_tilde2 = 110;
}
parameters {
  real<lower=0> epsilon_sigma;
  
  matrix<lower=0>[C, num_basis] alpha;
}
transformed parameters {
  matrix[C, T - 1] shock = rep_vector(0, C, T - 1);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  
  for (c in 1 : C) {
    for (t in 1 : (T - 1)) {
      transition_function[c, t] = rate_spline(y[c, t], P_tilde, P_tilde2,
                                              alpha[c,  : ], ext_knots,
                                              num_basis, spline_degree);
    }
  }
}
model {
  epsilon_sigma ~ std_normal();
  diff ~ normal(to_vector(transition_function) + to_vector(shock),
                epsilon_sigma);
  
  to_vector(alpha) ~ std_normal();
}
generated quantities {
  matrix[C, Tpred] eta;
  matrix[C, num_grid] transition_function_pred;
  matrix[C, Tpred] shock2 = rep_vector(C, Tpred);
  
  eta[1 : C, 1 : T] = y;
  
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;
  
  for (t in (T + 1) : Tpred) {
    for (c in 1 : C) {
      real transition = rate_spline(eta[c, t - 1], P_tilde, P_tilde2,
                                    alpha[c], ext_knots, num_basis,
                                    spline_degree);
      eta[c, t] = eta[c, t - 1] + transition
                  + error_rng(eta[c : c, t - 1], epsilon_params);
    }
  }
  
  for (c in 1 : C) {
    for (i in 1 : num_grid) {
      transition_function_pred[c, i] = rate_spline(grid[i], 0, 1, alpha[c],
                                                   ext_knots, num_basis,
                                                   spline_degree);
    }
  }
}

