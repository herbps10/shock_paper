cv_fit <- function(data, held_out, ...) {
	fpemplus(
		data,
		y = "contraceptive_use_modern",
		se = "se_modern",
		year = "year",
		source = "data_series_type",
		area = "name_country",
		held_out = held_out, 
		
		start_year = 1970, 
		end_year = 2030, 
		t_star = 1990,
		
		# Hierarchical setup
		hierarchical_level     = c("intercept", "name_region", "name_sub_region", "name_country"), 
		hierarchical_splines   = c("intercept", "name_region", "name_sub_region", "name_country"),
		hierarchical_asymptote = c("intercept", "name_country"),
		
		# Stan sampler settings
		parallel_chains = 4,
		refresh = 50, 
		...
	)
}

cv_fit_cutoff <- function(data, cutoff, ...) {
  held_out <- data$year >= cutoff
  cv_fit(data, held_out, ...) 
}
