library(tidyverse)
library(wpp2024)
#library(BayesTransitionModels)
library(tidybayes)
#library(bayesLife)

index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cache_path <- glue::glue("/gpfs/scratch/susmah01/shock_paper/validations-{index}.rds")

#if(file.exists(cache_path)) stop("File exists")

cmdstanr::set_cmdstan_path("/gpfs/data/diazi07lab/cmdstan-2.38.0/")

source("R/lifeplus.R")
source("R/process_lifeplus.R")

data(UNlocations, package = "wpp2024")
data(e0M1, package = "wpp2024")
data(pop1, package = "wpp2024")
data(include_2010, package = "bayesLife")

included_codes <- include_2010 |> filter(include_code %in% 1:2) |> pull(country_code)

large_countries <- pop1 |>
  filter(`2023` >= 1e3) |>
  pull(name)

datM <- e0M1 |>
  filter(name %in% large_countries) |>
  pivot_longer(cols = `1950`:`2023`, names_to = "period", values_to = "e0") |>
  filter(country_code %in% included_codes) |>
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2024")

#
# Validations
#
validation_cutoff <- function(model, cutoff_year, scale_global, outlier_threshold) {
  fit <- lifeplus(
    datM,
    y = "e0", 
    year = "year",
    area = "name",
    source = "source",
    start_year = 1950,
    end_year = max(datM$year),
    held_out = datM$year > cutoff_year,
    
    outlier_threshold = outlier_threshold,
    
    model = model,

    hierarchical = FALSE,
    centered = FALSE,

    adapt_delta = 0.999,
    max_treedepth = 15,
    parallel_chains = 8,
    chains = 8,
    iter_warmup = 500,
    iter_sampling = 500,

    output_dir = "/gpfs/scratch/susmah01/shock_paper/draws/",
    
    extra_stan_data = list(
      scale_global = scale_global,
      slab_scale = 10,
      slab_df = 6
    )
  )
  print(fit$samples$cmdstan_diagnose())
  return(fit)
}

validations <- expand_grid(
  cutoff_year = c(2013, 2018),
  model = c("logistic", "logistic_shock"),
  outlier_threshold = c(5, 1e3),
  scale_global = c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
) |>
  filter(!(model == "logistic_shock" & outlier_threshold == 5)) |>
  filter(!(model == "logistic" & scale_global != 1e-2))

print(nrow(validations))

validations <- validations[(1:nrow(validations) %% 18 + 1) == index, ]

print(glue::glue("Starting validation index {index} with {nrow(datM)} rows for {length(unique(datM$name))} countries"))
print(glue::glue("Cutoff: {validations$cutoff_year} model: {validations$model} outlier_threshold: {validations$outlier_threshold} scale_global: {validations$scale_global}"))

validations <- validations |>
  mutate(fit = pmap(list(model, cutoff_year, scale_global, outlier_threshold), validation_cutoff))

write_rds(validations, cache_path) 

