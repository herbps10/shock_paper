functions {
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }

  real shock_rng(real nu_local, real c_slab, real tau, int constrain_negative) {
    real shock_raw_pred;
    if(constrain_negative == 1) {
      shock_raw_pred = normal_lub_rng(0, 1, negative_infinity(), 0);
    }
    else {
      shock_raw_pred = normal_rng(0, 1);
    }
    real local_shrinkage_pred = student_t_rng(nu_local, 0, 1);
    real truncated_local_shrinkage_pred = sqrt(square(c_slab) * square(local_shrinkage_pred) ./ (square(c_slab) + square(tau) * square(local_shrinkage_pred)));
    return shock_raw_pred * truncated_local_shrinkage_pred * tau;
  }
}
data {
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  int<lower=0, upper=1> constrain_negative;
}
transformed data {
  shock_term = 1;
  generate_shock_free = 1;
  real tau = scale_global;
  int nu_local = 3;
  //int nu_global = 1;
}
parameters {
  vector<upper=(constrain_negative == 1 ? 0 : positive_infinity())>[C * T_shocks] shock_raw;
  vector<lower=0>[C * T_shocks] lambda;
  real<lower=0> caux;

  //real<lower=0> tau;
}
transformed parameters {
  real<lower=0> c_slab = slab_scale * sqrt(caux);
  vector<lower=0>[C * T_shocks] lambda_tilde = sqrt(c_slab^2 * square (lambda) ./ (c_slab^2 + tau^2 * square(lambda)));
  shock = to_matrix(shock_raw .* lambda_tilde * tau * epsilon_sigma, C, T_shocks);
}
model {
  shock_raw ~ std_normal();
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  lambda ~ student_t(nu_local, 0, 1);
  //tau ~ student_t(nu_global, 0, scale_global * epsilon_sigma);
}
generated quantities {
  real lambda_tilde_sd = sd(lambda_tilde);
  
  for(t in T:(T_shocks + Tpred - T)) {
    for(c in 1:C) {
      shock2[c, t] = shock_rng(nu_local, c_slab, tau, constrain_negative);
    }
  }
}
