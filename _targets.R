library(targets)
library(tarchetypes)
library(tidyverse)

source("R/fpem_cv.R")

tar_option_set(
  packages = c("tidyverse", "fpemplus")
)

library(future)
plan(multicore)

setup <- tribble(
	~model, ~scale_global, ~ssq
	#"spline", 0,
	"shock", 1e-2, 1,
	"shock", 1e-4, 1,
	"shock", 1e-6, 1,

	"shock", 1e-2, 10,
	"shock", 1e-4, 10,
	"shock", 1e-6, 10
)

fits_target <- tar_map(
	unlist = FALSE,
	values = setup,
	tar_target(fit, fpemplus(
    analysis_data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",
    
    start_year = 1970, 
    end_year = 2030, 
		t_star = 1990,
    
    model = ifelse(model == "spline", "spline", here::here("stan/fpem_spline_shock.stan")),
    
    smoothing = TRUE,
    
    # Spline setup
    spline_degree = 2,
    num_knots = 5,
    
    # Hierarchical setup
    hierarchical_level     = c("intercept", "name_region", "name_sub_region", "name_country"), 
    hierarchical_splines   = c("intercept", "name_region", "name_sub_region", "name_country"),
    hierarchical_asymptote = c("intercept", "name_country"),

    # Prior settings
    tau_prior = "normal(0, 2)",
    rho_prior = "uniform(0, 1)",

		extra_stan_data = list(
			scale_global = scale_global,
			ssq = ssq,
			slab_scale = 1,
			slab_df = 1
		),

		output_dir = "./output",
    
    # Stan sampler settings
    adapt_delta = 0.999,
    max_treedepth = 14,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 5,
    parallel_chains = 4,
    refresh = 50
  )),
	tar_target(fit_shock_cv_2010, cv_fit_cutoff(
    analysis_data,
		cutoff = 2010,
    
    model = ifelse(model == "spline", "spline", here::here("stan/fpem_spline_shock.stan")),
    
    smoothing = TRUE,
    
    # Spline setup
    spline_degree = 2,
    num_knots = 5,
    
    # Prior settings
    tau_prior = "normal(0, 2)",
    rho_prior = "uniform(0, 1)",

		extra_stan_data = list(
			scale_global = scale_global,
			ssq = ssq,
			slab_scale = 1,
			slab_df = 1
		),

		output_dir = "./output",
    
    # Stan sampler settings
    adapt_delta = 0.999,
    max_treedepth = 14,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 5
  )),
	tar_target(fit_shock_cv_2015, cv_fit_cutoff(
    analysis_data,
		cutoff = 2015,
    
    model = ifelse(model == "spline", "spline", here::here("stan/fpem_spline_shock.stan")),
    
    smoothing = TRUE,
    
    # Spline setup
    spline_degree = 2,
    num_knots = 5,
    
    # Prior settings
    tau_prior = "normal(0, 2)",
    rho_prior = "uniform(0, 1)",

		extra_stan_data = list(
			scale_global = scale_global,
			ssq = ssq,
			slab_scale = 1,
			slab_df = 1
		),

		output_dir = "./output",
    
    # Stan sampler settings
    adapt_delta = 0.999,
    max_treedepth = 14,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 5
  ))
)

list(
  tar_target(analysis_data, national_data(
    #countries = c("Somalia", "India", "Bangladesh", "Zimbabwe",
    #              "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala", 
		#							"Rwanda", "Timor-Leste", "Viet Nam", "Ghana", "Nigeria",
    #            	"Ethiopia", "Uganda"),
		countries = c(),
    start_year = 1970
  ) %>%
		filter(name_country != "Other non-specified areas") %>%
		filter(modern_method_bias == "None", has_geographical_region_bias == "N") %>%
		mutate(name_region = case_when(
			name_region == "Europe" ~ "Europe and Northern America",
			name_region == "Northern America" ~ "Europe and Northern America",
			TRUE ~ name_region
		))
	),
	fits_target
)

