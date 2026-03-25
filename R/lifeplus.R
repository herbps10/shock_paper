library(splines)

#' @import tibble
hierarchical_model_matrix <- function(columns, data) {
  mat <- matrix(NA, nrow = nrow(data), ncol = 0)
  assign <- c()
  
  index <- tibble::tibble(i = numeric(0), column = character(0), level = character(0))
  
  i <- 1
  if("intercept" %in% columns) {
    mat <- cbind(mat, rep(1, nrow(mat)))
    assign <- c(assign, 0)
    
    index[i, ] <- tibble(i = i, column = "intercept", level = "intercept")
    
    i <- i + 1
  }
  
  for(column in columns) {
    for(l in levels(factor(data[[column]]))) {
      mat <- cbind(mat, as.numeric(data[[column]] == l))
      assign <- c(assign, which(column == columns))
      index[i, ] <- tibble::tibble(i = i, column = column, level = l)
      
      i <- i + 1
    }
  }
  
  list(
    assign = assign,
    matrix = mat,
    index = index
  )
}

#' @import purrr
hierarchical_data <- function(data, hierarchy) {
  model_matrix <- hierarchical_model_matrix(hierarchy, data)
  n_terms <- ncol(model_matrix$mat)
  re <- unique(model_matrix$assign)
  n_re <- length(re)
  re_start <- map_int(re, function(x) min(which(model_matrix$assign == x)))
  re_end   <- map_int(re, function(x) max(which(model_matrix$assign == x)))
  
  list(
    model_matrix = model_matrix,
    n_terms = n_terms,
    n_re = n_re,
    re_start = re_start,
    re_end = re_end
  )
}

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
  transition = "logistic",
  shock = FALSE,
  data_model = "normal",
  
  num_knots = 7,
  spline_degree = 2,
  outlier_threshold = 1000,
  
  hierarchical = TRUE,
  centered = TRUE,
  
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
  root <- rprojroot::is_git_root                                                                         
  basepath <- root$find_file("stan")  
  
  stan_file_path <- paste0(
    basepath, "/",
    paste0(c(transition, data_model, ifelse(shock == TRUE, "shock", "noshock")), collapse = "_"),
    ".stan"
  )
  
  stan_model <- cmdstanr::cmdstan_model(
    stan_file_path,
    dir = tempdir()
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
  
  a_data       <- hierarchical_data(country_index, hierarchical_splines)
  
  # Set up spline basis
  knots <- sort(c(seq(0, max(data[[y]]) / 110, length.out = num_knots), 1, 2))
  
  grid <- c(seq(from = 0, to = 110, by = 1)) # generating inputs
  num_grid <- length(grid)
  
  if(length(held_out) == 1 && held_out == FALSE) {
    obs <- data |> select(t, c, e0) |> pivot_wider(names_from = "t", values_from = "e0") |> select(-c) |> as.matrix()
  }
  else {
    obs <- data[held_out == 0,] |> select(t, c, e0) |> pivot_wider(names_from = "t", values_from = "e0") |> select(-c) |> as.matrix()
  }
  
  B <- t(bs(grid, knots = knots, degree = spline_degree, intercept = FALSE))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  ext_knots <- c(rep(knots[1], spline_degree), knots, rep(knots[length(knots)], spline_degree))
  
  a_lower_bound <- 0.01
  a_upper_bound <- 10 
  
  stan_data <- c(extra_stan_data, list(
    C = nrow(obs),
    T = ncol(obs),
    Tpred = max(time_index$t),
    
    y = obs,
    
    hierarchical = as.numeric(hierarchical),
    centered = as.numeric(centered),
    
    outlier_threshold = outlier_threshold,
    
    num_grid = num_grid,
    grid = grid,
    
    # Spline settings
    num_knots = length(knots),
    knots = knots,
    
    spline_degree = spline_degree,
    B = B,
    
    Delta1_constrain = 1, Delta1_lower = 0, Delta1_upper = 50,  Delta1_prior_mean = 0, Delta1_prior_sd = 1,
    Delta2_constrain = 1, Delta2_lower = 0, Delta2_upper = 50,  Delta2_prior_mean = 0, Delta2_prior_sd = 1,
    Delta3_constrain = 1, Delta3_lower = 0, Delta3_upper = 50,  Delta3_prior_mean = 0, Delta3_prior_sd = 1,
    Delta4_constrain = 1, Delta4_lower = 5, Delta4_upper = 50,  Delta4_prior_mean = 0, Delta4_prior_sd = 1,
    k_constrain = 1,      k_lower = 0,      k_upper = 10,       k_prior_mean = 0,      k_prior_sd = 1,
    z_constrain = 1,      z_lower = 0,      z_upper = 1.15/5,   z_prior_mean = 0,      z_prior_sd = 1,
    
    alpha_constrain = 1,  alpha_lower = 0,  alpha_upper = 10,   alpha_prior_mean = -2, alpha_prior_sd = 2,
    beta_constrain = 0,   beta_lower = 0,   beta_upper = 1,     beta_prior_mean = 0,   beta_prior_sd = 1
  ))
    
  start <- Sys.time()
  fit <- stan_model$sample(
    stan_data,
    save_latent_dynamics = TRUE,
    ...
  )
  elapsed <- Sys.time() - start
  
  result <- list(samples = fit,
                 data = original_data,
                 stan_data = stan_data,
                 time_index = time_index,
                 country_index = country_index,
                 
                 elapsed = elapsed,
                 
                 # Save arguments
                 y = y,
                 year = year,
                 source = source,
                 area = area,
                 held_out = held_out,
                 
                 transition = transition,
                 shock = shock,
                 data_model = data_model)
  
  cat("Extracting posteriors...\n")
  
  result$posteriors <- process_life_fit(result, ifelse(is.null(args$parallel_chains), 1, args$parallel_chains))


  result$diagnose <- fit$diagnostic_summary()
  
  attr(result, "class") <- "fpemplus"
  
  result
}

