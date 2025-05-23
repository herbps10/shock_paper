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
#' @param hierarchical_splines vector specifying hierarchical structure for spline coefficients (see Details).
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
  
  country_specific_global_shrinkage = FALSE,
  
  extra_stan_data = list(),
  
  hierarchical_splines = c("intercept", area),
  
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
  
  if(length(hierarchical_splines) == 0) {
    stop("No hierarchical structure supplied for the spline coefficients. See the hierarchical_splines argument.")
  }
  
  # Make sure there are no NAs in supplied columns
  BayesTransitionModels:::check_nas(data, y)
  BayesTransitionModels:::check_nas(data, year)
  
  # Initialize start and end year if necessary
  if(is.na(start_year)) start_year <- min(data[[year]])
  if(is.na(end_year)) end_year <- max(data[[year]])
  
  # Make sure the observed data are within the estimation period
  if(sum(!(data[[year]] %in% start_year:end_year)) > 0) {
    stop(glue::glue("Observations included in dataset that fall outside the estimation period ({start_year} to {end_year})."))
  }
  
  ###### Load model #####
  include_paths <- system.file("include", package = "BayesTransitionModels")
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
  else {
    stop(glue::glue("Model {model} not supported. Currently \"spline\" is the only supported model."))
  }
  
  stan_model <- cmdstanr::cmdstan_model(
    stan_file_path,
    include_paths = include_paths
  )
  
  #
  # Setup data for Stan
  #
  
  # Create district index for matching district and district index
  hierarchical_column_names <- unique(c(
    hierarchical_splines
  )) %>%
    setdiff("intercept")
  
  # Make sure there are no NAs in any of the columns
  for(column in hierarchical_column_names) {
    if(column == "intercept") next
    BayesTransitionModels:::check_nas(data, column)
  }
  
  country_index <- data %>%
    dplyr::distinct(!!! syms(hierarchical_column_names)) %>%
    dplyr::mutate(c = 1:n())
  
  # Create year lookup table
  time_index <- tibble(
    year = seq(start_year, end_year, 5),
    t = 1:length(year)
  ) 
  
  year_by <- c()
  year_by[year] = year
  data <- data %>%
    dplyr::left_join(time_index, by = year_by) %>%
    dplyr::left_join(country_index, by = hierarchical_column_names)
  
  if(length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  }
  else {
    held_out = as.numeric(held_out)
  }
  
  t_last <- max(data$t)
  
  # Set up hierarchical structures
  a_data       <- BayesTransitionModels:::hierarchical_data(country_index, hierarchical_splines)
  
  # Set up spline basis
  knots <- sort(c(seq(0, 1, length.out = num_knots), 1000))
  grid <- c(seq(from = 0, to = 1, by = .05), 1000) # generating inputs
  
  B <- t(bs(grid, knots = knots, degree = spline_degree, intercept = FALSE))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  ext_knots <- c(rep(knots[1], spline_degree), knots, rep(knots[length(knots)], spline_degree))
  
  a_lower_bound <- 0.01
  a_upper_bound <- 10 
  
  stan_data <- c(extra_stan_data, list(
    C = nrow(country_index),
    T = nrow(time_index),
    N = nrow(data),
    held_out = held_out,
    t_last = t_last,
    
    time = array(data$t),
    country = array(data$c),
    
    y = array(data[[y]]),
    
    a_n_terms = a_data$n_terms,
    a_n_re = a_data$n_re,
    a_re_start = array(a_data$re_start),
    a_re_end = array(a_data$re_end),
    a_model_matrix = a_data$model_matrix$mat,
    
    # Spline settings
    num_knots = length(knots),
    knots = knots,
    
    num_grid = num_grid,
    spline_degree = spline_degree,
    grid = grid,
    B = B,
    
    a_lower_bound = a_lower_bound,
    a_upper_bound = a_upper_bound,
    R = R,
    
    crisis_in_projections = crisis_projections,
    country_specific_global_shrinkage = country_specific_global_shrinkage
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
                 hierarchical_splines = hierarchical_splines,
                 held_out = held_out,
                 model = model)
  
  cat("Extracting posteriors...\n")
  
  result$posteriors <- process_life_fit(result, ifelse(is.null(args$parallel_chains), 1, args$parallel_chains))
  
  attr(result, "class") <- "fpemplus"
  
  result
}

