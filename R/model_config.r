library(jsonlite)

generate_config <- function(transition_model, data_model, shocks) {
  if(!(transition_model %in% c("logistic", "spline", "gp"))) stop("transition_model must be one of: logistic, spline, gp")
  if(!(data_model %in% c("normal", "outlier", "mixture"))) stop("data_model must be one of: normal, outlier, mixture")
  if(!(shocks %in% c(FALSE, TRUE))) stop("shocks must be one of: FALSE, TRUE")
  
  path <- glue::glue("../{transition_model}_{data_model}_{ifelse(shocks, 'shock', 'noshock')}.stan")
  
  deps <- list(
    output = path, 
    "base.stan" = c()
  )
  
  if(transition_model == "gp" || data_model == "heteroskedastic") {
    deps <- c(deps, list(
      "modules/approximate_gp.stan" = c()
    ))
  }
  
  if(shocks == TRUE) {
    deps <- c(deps, list(
      "modules/shock.stan" = c()
    ))
  }
  
  if(data_model == "normal") {
    deps <- c(deps, list(
      "modules/data_model_normal.stan" = c()
    ))
  }
  else if(data_model == "outlier") {
    deps <- c(deps, list(
      "modules/data_model_outlier.stan" = c()
    ))
  }
  else if(data_model == "mixture") {
    deps <- c(deps, list(
      "modules/data_model_mixture.stan" = c()
    ))
  }
  else if(data_model == "heteroskedastic") {
    deps <- c(deps, list(
      "modules/data_model_heteroskedastic.stan" = c()
    ))
  }
  
  if(transition_model == "logistic") {
    deps <- c(deps, list(
      "modules/Delta.stan" = c(),
      "modules/hierarchical_matrix_cholesky.stan" = list("var" = "Delta", num = "D"),
      "modules/transition_double_logistic.stan" = c()
    ))
  }
  else if(transition_model == "gp") {
    deps <- c(deps, list(
      "modules/hierarchical_matrix.stan" = list("var" = "beta", "num" = "M"),
      "modules/transition_gp.stan" = c()
    ))
  }
  else if(transition_model == "spline") {
    deps <- c(deps, list(
      "modules/hierarchical_matrix.stan" = list("var" = "alpha", "num" = "num_basis"),
      "modules/transition_spline.stan" = c()
    ))
  }
  
  
  
  deps
}

configs <- expand_grid(
  transition_model = c("logistic", "spline", "gp"),
  data_model = c("normal", "outlier", "mixture"),
  shock = c(FALSE, TRUE)
) |>
  mutate(config = pmap(list(transition_model, data_model, shock), generate_config))
  
jsonlite::toJSON(configs$config, auto_unbox = TRUE) |>
  str_replace_all("\\.[0-9]", "") |>
  write_file("stan/config/config.json")
