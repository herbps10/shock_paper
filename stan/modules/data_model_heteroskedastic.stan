functions {
  real error_rng(vector x, vector params, int shock) {
    int M = size(params) - 1;
    real L = params[M + 1];
    matrix[1, M] p;
    for(m in 1:M) p[, m] = phi(L, m, x / 110);

    return (exp(p * params[1:M]))[1];
  }
}
parameters {
  real<lower=0> epsilon_rho;
  real<lower=0> epsilon_scale;
  vector[M]     epsilon_beta;
}
transformed parameters {
  vector[M] epsilon_diagSPD;
  vector[M] epsilon_SPD_beta;
  vector[C * (T - 1)] epsilon_sigma;

  for(m in 1:M) epsilon_diagSPD[m] = sqrt(spd(epsilon_scale, epsilon_rho, sqrt(flambda(L, m))));
  epsilon_SPD_beta = epsilon_diagSPD .* epsilon_beta;
  epsilon_sigma = exp(PHI * epsilon_SPD_beta);
}
model {
  epsilon_scale ~ std_normal();
  epsilon_rho   ~ inv_gamma(5, 5);
  epsilon_beta  ~ std_normal();

  diff ~ normal(to_vector(transition_function) + to_vector(shock), epsilon_sigma);
}
generated quantities {
  vector[M + 1] epsilon_params;
  epsilon_params[1:M] = epsilon_SPD_beta;
  epsilon_params[M + 1] = L;

  vector[num_grid] epsilon_sigma_pred = exp(PHI_grid * to_vector(epsilon_SPD_beta));
}
