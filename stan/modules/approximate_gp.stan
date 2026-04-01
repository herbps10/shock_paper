functions {
  real flambda(real L, int m) {
		real lam;
		lam = ((m*pi())/(2*L))^2;
				
		return lam;
	}
	real spd(real alpha, real rho, real w) {
		return (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2).*(w^2));
	}
	
	vector phi(real L, int m, vector x) {
		vector[rows(x)] fi;
		fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
				
		return fi;
	}
}
data {
  real L;
  int<lower=1> M;
}
transformed data {
  matrix[C * (T - 1), M] PHI;
  matrix[num_grid, M] PHI_grid;
  for(m in 1:M) {
    PHI_grid[, m] = phi(L, m, grid / 110);
    for(c in 1:C) {
      int start = (c - 1) * (T - 1) + 1;
      int end = c * (T - 1);
      PHI[start:end, m] = phi(L, m, to_vector((y[c, 1:(T - 1)]) / 110));
    }
  }
}
