# -------------------------------
# Functions
# -------------------------------

load_libraries <- function() {
  if (!require("tidyverse")) install.packages("tidyverse"); require ("tidyverse")
  if (!require("foreach")) install.packages("foreach"); require ("foreach")
  if (!require("ggplot2")) install.packages("ggplot2"); require ("ggplot2")
  if (!require("bayesplot")) install.packages("bayesplot"); require ("bayesplot")
  if (!require("brms")) install.packages("brms"); require ("brms")
  if (!require("rstan")) install.packages("rstan"); require ("rstan")
  if (!require("caret")) install.packages("caret"); require ("caret")
  if (!require("rsample")) install.packages("rsample"); require ("rsample")
  if (!require("tidybayes")) install.packages("tidybayes"); require ("tidybayes")
  if (!require("loo")) install.packages("loo"); require ("loo")
  if (!require("furrr")) install.packages("furrr"); require ("furrr")
  if (!require("lme4")) install.packages("lme4"); require ("lme4")
}

load_posthoc_libraries <- function() {
  if (!require("lmerTest")) install.packages("lmerTest"); require ("lmerTest")
  if (!require("kableExtra")) install.packages("kableExtra"); require ("kableExtra")
  if (!require("scales")) install.packages("scales"); require ("scales")
}

cent_scale <- function(dat, method = c("GM", "CWC")) {
  
  method <- match.arg(method)
  
  if (method == "GM") {
    
    # grand mean centering
    df <- dat %>%
      mutate_if(is.numeric, ~((.x - mean(.x, na.rm=TRUE)) / sd(.x, na.rm=TRUE)))
    
  } else if (method == "CWC") {
    
    # scale numeric data within clusters
    df <- dat %>%
      group_by(user_id) %>%
      mutate_if(is.numeric, ~((.x - mean(.x, na.rm=TRUE)) / sd(.x, na.rm=TRUE))) %>%
      mutate_if(is.numeric, ~if_else(is.nan(.), 0, .)) %>%
      ungroup()
    
  } else {
    print("ERROR")
  }
  
  return(df)
}

modeling_func <- function(train, key, chain_num = 4, cores = 4,
                          iter_num = 2000, mod = c("all", "WP", "int")) {
  
  mod <- match.arg(mod)
  
  # create formula to evaluate all variables in relation to cog outcomes
  df_all <- train %>% select(-value, -user_id, -study_time, -sitting_id)
  formula_all <- paste0(names(df_all), collapse=" + ")
  
  # create formula to evaluate WP predictors of variation around WP means
  to_elim <- c("age", "education", "gender", "englishPrimary",
               "wakeUpTypical", "wakeUpEarliest", "wakeUpLatest",
               "goToSleepTypical", "goToSleepEarliest", "goToSleepLatest",
               "user_id", "value", "study_time", "screenW", "screenH", "sitting_id")
  df_WP <- train %>% select(-all_of(to_elim))
  formula_WP <- paste0(names(df_WP), collapse=" + ")
  
  # if - else logic
  if (mod == "all") {
    formula <- as.formula(paste0("value ~ 1 + ", formula_all, 
                                 " + age:cos_studytime + age:sin_studytime + I(age^2) + (1 | user_id)"))
    fit_bayes_horse <- brm(formula, 
                           data = train, family = 'gaussian',
                           prior = c(prior(horseshoe(df = 4), class = "b"),
                                     prior(cauchy(0, 1), class = "sigma"),
                                     prior(normal(0, 1), class = "sd"),
                                     prior(normal(0, 1), class = "Intercept")),
                           chains = chain_num, iter = iter_num,
                           control = list(adapt_delta = .99, max_treedepth = 15),
                           cores = 4, file = paste0(path, key, "_re"))
  } else if (mod == "WP") {
    formula <- as.formula(paste0("value ~ 1 + ", formula_WP, " + (1 | user_id)"))
    fit_bayes_horse <- brm(formula, 
                           data = train, family = 'gaussian',
                           prior = c(prior(horseshoe(df = 4), class = "b"),
                                     prior(cauchy(0, 1), class = "sigma"),
                                     prior(normal(0, 1), class = "sd"),
                                     prior(normal(0, 1), class = "Intercept")),
                           chains = chain_num, iter = iter_num,
                           control = list(adapt_delta = .99, max_treedepth = 15),
                           cores = 4, file = paste0(path, key, "_fe"))
  } else if (mod == "int") {
    formula <- as.formula(paste0("value ~ 1 + (1 | user_id)"))
    fit_bayes_horse <- brm(formula, 
                           data = train, family = 'gaussian',
                           prior = c(prior(cauchy(0, 1), class = "sigma"),
                                     prior(normal(0, 1), class = "sd"),
                                     prior(normal(0, 1), class = "Intercept")),
                           chains = chain_num, iter = iter_num,
                           control = list(adapt_delta = .99, max_treedepth = 15), 
                           cores = 4, file = paste0(path, key, "_int"))
  } else {
    print("ERROR")
  }
  
  return(fit_bayes_horse)
}

rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  return(sqrt(mean((yrep_mean - y)^2)))
}

R2 <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  x <- cor.test(y, yrep_mean)
  corre <- x$estimate**2
  return(corre[[1]])
}
