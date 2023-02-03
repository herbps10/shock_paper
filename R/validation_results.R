library(tidyverse)
library(targets)
library(tidybayes)

#tar_load(starts_with("fit_shock_cv"))
tar_load(fit_shock_cv_2010_shock_0.01)
tar_load(fit_shock_cv_2010_shock_1e.04)
tar_load(fit_shock_cv_2010_shock_1e.06)
tar_load(fit_shock_cv_2010_spline_0)

validation_measures <- function(fit) {
	fit$posteriors$generated_quantities %>% 
		median_qi(y_pred) %>% 
		bind_cols(fit$data) %>% 
		mutate(
			below = contraceptive_use_modern < .lower,
			above = contraceptive_use_modern > .upper,
			covered = .lower <= contraceptive_use_modern, .upper >= contraceptive_use_modern, 
			ci_width = .upper - .lower,
			error = contraceptive_use_modern - y_pred,
			held_out = fit$held_out
		) %>%
		filter(held_out == 1) %>%
		group_by(name_country) %>%
		filter(year == max(year)) %>%
		ungroup() %>%
		summarize(
			below = mean(below),
			coverage = mean(covered), 
			above = mean(above),
			ci_width = median(ci_width),
			median_error = median(error),
			median_abs_error = median(abs(error))
		)
}

results <- bind_rows(
	validation_measures(fit_shock_cv_2010_shock_0.01) %>% mutate(model = "with shocks 1e-2"),
	validation_measures(fit_shock_cv_2010_shock_1e.04) %>% mutate(model = "with shocks 1e-4"),
	validation_measures(fit_shock_cv_2010_shock_1e.06) %>% mutate(model = "with shocks 1e-6"),
	validation_measures(fit_shock_cv_2010_spline_0) %>% mutate(model = "without shocks")
)
