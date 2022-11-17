library(invgamma)

simulate_horseshoe <- function(N, scale_global = 1) {
  tau <- rt(N, df = 1) * scale_global
  lambda <- rcauchy(N)
  
  beta <- rnorm(N, mean = 0, sd = abs(tau * lambda))
  
  beta
}

simulate_regularized_horseshoe <- function(N, scale_global = 1, local_df, slab_df = 1) {
  tau <- rt(1, df = 1) * scale_global
  lambda <- rt(N, df = local_df)
  csq <- rinvgamma(1, 0.5 * slab_df, 0.5 * slab_df)
  
  lambda_tilde = csq * lambda^2 / (csq + tau^2 * lambda^2)
  beta <- rnorm(N, 0, abs(tau) * sqrt(lambda_tilde))
  beta
}
threshold <- 0.1
results <- rerun(1e3,
  simulate_regularized_horseshoe(8e5, scale_global = 0.0001, slab_df = 1, local_df = 1)
)

p <- mean(map_dbl(results, function(x) mean(x > threshold)))

(1 - p) * 100

mean(unlist(results))
hist(unlist(results))

x <- 
hist(x)

mean(x > 1)
hist(x[abs(x) < 1])

