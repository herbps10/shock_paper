process_life_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1
  
  #
  # eta and epsilon summaries
  #
  temporal_variables <- c("eta")
  
  temporal <- fit$samples$summary(temporal_variables, ~stats::quantile(.x, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)), .cores = parallel_chains) %>%
    tidyr::separate(.data$variable, c("variable", "index"), "\\[") %>%
    dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) %>%
    tidyr::separate(.data$index, c("c", "t"), ",") %>%
    dplyr::mutate_at(vars(c, t), as.integer) %>%
    dplyr::left_join(fit$country_index, by = "c") %>%
    dplyr::left_join(fit$time_index, by = "t")
  
  #
  # Transition function summaries
  #
  #transition <- list()
  #transition_quantiles <- list()
  #for(column in fit$hierarchical_splines) {
  #  transition[[column]] <- BayesTransitionModels:::extract_rate_vs_level_subhierarchical(fit, fit$hierarchical_splines, column, "after")
  #  transition_quantiles[[column]] <- transition[[column]] %>%
  #    dplyr::group_by(.data$name, .data$x) %>%
  #    tidybayes::median_qi(.data$Y, .width = c(0.5, 0.8, 0.95))
  #}
  
  
  #
  # Hierarchical distributions
  #
  a_sigma <- fit$samples$draws(c("a_sigma")) %>%
    tidybayes::spread_draws(a_sigma[i, j])
  
  ans <- list(
    temporal = temporal,
    #transition = transition,
    #transition_quantiles = transition_quantiles,
    a_sigma = a_sigma
  )
  
  ans
}
