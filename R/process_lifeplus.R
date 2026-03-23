process_life_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1
  
  #
  # eta and epsilon summaries
  #
  temporal_variables <- c("eta")
  if(fit$shock == TRUE) temporal_variables <- c(temporal_variables, "shock2", "eta_shockfree")
  if(fit$data_model == "mixture") temporal_variables <- unique(c(temporal_variables, "eta_shockfree"))
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

  transition_functions <- fit$samples$draws("transition_function_pred") |> spread_draws(transition_function_pred[c, i]) |>
    left_join(fit$country_index) |>
    left_join(tibble(i = 1:length(fit$stan_data$grid), x = fit$stan_data$grid)) |>
    filter(x < 1000) |>
    group_by(name, x) |>
    median_qi(transition_function_pred, .width = c(0.5, 0.8, 0.95))
  
  transition_function_mean <- NULL
  if(fit$stan_data$hierarchical == 1) {
    transition_function_mean <- fit$samples$draws("transition_function_pred_mean") |> spread_draws(transition_function_pred_mean[i]) |>
      left_join(tibble(i = 1:length(fit$stan_data$grid), x = fit$stan_data$grid)) |>
      filter(x < 1000) |>
      group_by(x) |>
      median_qi(transition_function_pred_mean, .width = c(0.5, 0.8, 0.95))
  }
  
  transition_params <- NULL
  if(fit$transition == "logistic") {
    transition_params <- fit$samples$summary(c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z"), ~stats::quantile(.x, probs = c(0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999)), .cores = parallel_chains)  |>
      mutate_at(vars(ends_with("%")), as.numeric) |>
      tidyr::separate(.data$variable, c("variable", "index"), "\\[") |>
      dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) |>
      tidyr::separate(.data$index, c("c"), ",") |>
      dplyr::mutate_at(vars(c), as.integer) |>
      dplyr::left_join(fit$country_index, by = "c")
  }
  
  ans <- list(
    temporal = temporal,
    transition_params = transition_params,
    transition_functions = transition_functions,
    transition_function_mean = transition_function_mean
  )
  
  ans
}
