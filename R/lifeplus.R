library(splines)

#' Fit Lifeplus model
#'
#' @param data a data frame.
#' @param start_year start year of estimates.
#' @param end_year end year of estimates.
#' @param y column name of outcome.
#' @param year column name of outcome year.
#' @param area column name of the area of each observation
#' @param t_star reference year.
#' @param model which model to fit. Currently only "spline" is supported.
#' @param num_knots number of spline knots.
#' @param spline_degree spline degree. Degree 2 or 3 is supported.
#' @param held_out binary vector indicating which observations are held out. Set to FALSE to hold out no observations.
#' @param ... additional arguments for CmdStanModel::sample.
#'
#' @export
lifeplus <- function(
    # Data and column names
  data,
  y,
  year,
  area,
  source,
  
  start_year = NA,
  end_year = NA,
  t_star = NULL,
  
  crisis_projections = TRUE,
  
  # Model settings
  model = "spline",
  num_knots = 7,
  spline_degree = 2,
  outlier_threshold = 1000,
  
  hierarchical = TRUE,
  
  country_specific_global_shrinkage = FALSE,
  
  extra_stan_data = list(),
  
  # Out-of-sample validation
  held_out = FALSE,
  
  R = 1e3,
  
  # Stan settings
  ...
) {
  args <- list(...)
  
  # Save original dataset
  original_data <- data
  
  ###### Initial argument checks #####
  stopifnot(is.numeric(spline_degree))
  stopifnot(is.numeric(num_knots))
  
  if(!(spline_degree %in% c(2, 3))) {
    stop("spline_degree must be either 2 or 3.")
  }
  
  if(num_knots <= 0) {
    stop("num_knots must be greater than zero.")
  }
  
  if(nrow(data) == 0) {
    stop("Data has no rows.")
  }
  
  if(length(held_out) > 1) {
    if(length(held_out) != nrow(data)) stop(glue::glue("held_out (length {length(held_out)}) must be same size as dataset ({nrow(data)} rows)."))
  }
  
  if(start_year > end_year) {
    stop("start_year must be less than end year")
  }
  
  # Make sure there are no NAs in supplied columns
  #BayesTransitionModels:::check_nas(data, y)
  #BayesTransitionModels:::check_nas(data, year)
  
  # Initialize start and end year if necessary
  if(is.na(start_year)) start_year <- min(data[[year]])
  if(is.na(end_year)) end_year <- max(data[[year]])
  
  # Make sure the observed data are within the estimation period
  if(sum(!(data[[year]] %in% start_year:end_year)) > 0) {
    stop(glue::glue("Observations included in dataset that fall outside the estimation period ({start_year} to {end_year})."))
  }
  
  ###### Load model #####
  #include_paths <- system.file("include", package = "BayesTransitionModels")
  #stan_file_path <- system.file("stan/tfr_spline.stan", package = "BayesTransitionModels")
  
  if(model == "spline") {
    stan_file_path <- "stan/life_spline.stan"
  }
  else if(model == "shock") {
    stan_file_path <- "stan/life_spline_shock.stan"
  }
  else if(model == "shock2") {
    stan_file_path <- "stan/life_spline_shock2.stan"
  }
  else if(model == "logistic") {
    stan_file_path <- "stan/life_double_logistic.stan"
  }
  else if(model == "logistic_shock") {
    stan_file_path <- "stan/life_double_logistic_shock.stan"
  }
  else {
    stop(glue::glue("Model {model} not supported. Currently \"spline\" is the only supported model."))
  }
  
  stan_model <- cmdstanr::cmdstan_model(
    stan_file_path,
    dir = tempdir()
    #include_paths = include_paths
  )
  
  #
  # Setup data for Stan
  #
  
  country_index <- data |>
    dplyr::distinct(!!! syms(area)) |>
    dplyr::mutate(c = 1:n())
  
  # Create year lookup table
  time_index <- tibble(
    year = seq(start_year, end_year, 1),
    t = 1:length(year)
  ) 
  
  year_by <- c()
  year_by[year] = year
  data <- data |>
    dplyr::left_join(time_index, by = year_by) |>
    dplyr::left_join(country_index, by = area)
  
  if(length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  }
  else {
    held_out = as.numeric(held_out)
  }
  
  t_last <- max(data$t)
  
  # Set up spline basis
  if(model == "logistic" || model == "logistic_shock") {
    grid <- c(seq(from = 0, to = 110, by = 5)) # generating inputs
  }
  else {
    grid <- c(seq(from = 0, to = 1, by = .05)) # generating inputs
  }
  num_grid <- length(grid)
  
  if(length(held_out) == 1 && held_out == FALSE) {
    obs <- data |> select(t, c, e0) |> pivot_wider(names_from = "t", values_from = "e0") |> select(-c) |> as.matrix()
  }
  else {
    obs <- data[held_out == 0,] |> select(t, c, e0) |> pivot_wider(names_from = "t", values_from = "e0") |> select(-c) |> as.matrix()
  }
  
  stan_data <- c(extra_stan_data, list(
    C = nrow(obs),
    T = ncol(obs),
    Tpred = max(time_index$t),
    
    y = obs,
    
    hierarchical = as.numeric(hierarchical),
    
    outlier_threshold = outlier_threshold,
    
    num_grid = num_grid,
    grid = grid
  ))
  
  fit <- stan_model$sample(
    stan_data,
    save_latent_dynamics = TRUE,
    ...
  )
  
  result <- list(samples = fit,
                 data = original_data,
                 stan_data = stan_data,
                 time_index = time_index,
                 country_index = country_index,
                 
                 # Save arguments
                 y = y,
                 year = year,
                 source = source,
                 area = area,
                 held_out = held_out,
                 model = model)
  
  cat("Extracting posteriors...\n")
  
  result$posteriors <- process_life_fit(result, ifelse(is.null(args$parallel_chains), 1, args$parallel_chains))
  
  attr(result, "class") <- "fpemplus"
  
  result
}

