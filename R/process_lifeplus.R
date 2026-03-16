process_life_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1
  
  #
  # eta and epsilon summaries
  #
  if(fit$model == "spline" || fit$model == "logistic" || fit$model == "gp") {
    temporal_variables <- c("eta")
  } 
  else if(fit$model == "logistic_mixture") {
    temporal_variables <- c("eta", "eta_crisisfree")
  }
  else if(fit$model == "logistic_shock" || fit$model == "shock2") {
    temporal_variables <- c("eta", "eta_crisisfree", "shock2")
  }
  
  temporal <- fit$samples$summary(temporal_variables, ~stats::quantile(.x, probs = c(0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999)), .cores = parallel_chains) |>
    mutate_at(vars(ends_with("%")), as.numeric) |>
    tidyr::separate(.data$variable, c("variable", "index"), "\\[") |>
    dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) |>
    tidyr::separate(.data$index, c("c", "t"), ",") |>
    dplyr::mutate_at(vars(c, t), as.integer) |>
    dplyr::left_join(fit$country_index, by = "c") |>
    dplyr::left_join(fit$time_index, by = "t")
  
  #
  # Transition function summaries
  #
  #transition_function_mean <- fit$samples$draws("transition_function_mean") |> spread_draws(transition_function_mean[i]) |>
  #  left_join(tibble(i = 1:length(fit$stan_data$grid), x = fit$stan_data$grid)) |>
  #  filter(x < 1000) |>
  #  group_by(x) |>
  #  median_qi(transition_function_mean, .width = c(0.5, 0.8, 0.95))

  transition_functions <- fit$samples$draws("transition_function_pred") |> spread_draws(transition_function_pred[c, i]) |>
    left_join(fit$country_index) |>
    left_join(tibble(i = 1:length(fit$stan_data$grid), x = fit$stan_data$grid)) |>
    filter(x < 1000) |>
    group_by(name, x) |>
    median_qi(transition_function_pred, .width = c(0.5, 0.8, 0.95))
  
  #
  # Hierarchical distributions
  #
  #a_sigma <- fit$samples$draws(c("a_sigma")) |>
  #  tidybayes::spread_draws(a_sigma[i])
  
  ans <- list(
    temporal = temporal,
    transition_functions = transition_functions
    #transition_function_mean = transition_function_mean
    #a_sigma = a_sigma
  )
  
  ans
}
