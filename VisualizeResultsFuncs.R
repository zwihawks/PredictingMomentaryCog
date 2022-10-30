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
  if (!require("ggpubr")) install.packages("ggpubr"); require ("ggpubr")
  if (!require("ggrepel")) install.packages("ggrepel"); require ("ggrepel")
  if (!require("ggforce")) install.packages("ggforce"); require ("ggforce")
  if (!require("ggpmisc")) install.packages("ggpmisc"); require ("ggpmisc")
  if (!require("ggtext")) install.packages("ggtext"); require ("ggtext")
  if (!require("kableExtra")) install.packages("kableExtra"); require ("kableExtra")
  if (!require("bayestestR")) install.packages("bayestestR"); require ("bayestestR")
  if (!require("cowplot")) install.packages("cowplot"); require ("cowplot")
  if (!require("corrr")) install.packages("corrr"); require ("corrr")
  if (!require("ggcorrplot")) install.packages("ggcorrplot"); require ("ggcorrplot")
  if (!require("corrplot")) install.packages("corrplot"); require ("corrplot")
  if (!require("performance")) install.packages("performance"); require ("performance")
}

rmse <- function(y, yrep) {
  return(sqrt(mean((yrep - y)^2)))
}

R2 <- function(y, yrep) {
  x <- cor.test(y, yrep)
  cor <- x$estimate[[1]]**2
  lower <- x$conf.int[[1]]**2
  upper <- x$conf.int[[2]]**2
  return(c(cor, lower, upper))
}

unnest_kfold <- function(kfld, folds, widths) {
  for (i in 1:folds) {
    tmp <- kfld$fits[i,]$fit %>%
      gather_draws(., `b_.*`, regex = TRUE) %>% 
      mean_qi(.width = widths) %>% 
      mutate(fold = i)
    
    if (i == 1) {
      df <- tmp
    } else {
      df <- df %>% bind_rows(tmp)
    }
  }
  return(df)
}

perf_cal <- function(kf, kfp, method = c("avg", "ind"), metric = c("R2", "RMSE")) {
  method <- match.arg(method)
  
  pr <- kf$fits[, "predicted"]
  n_folds <- length(pr)
  post_draws <- nrow(kfp$yrep)
  
  if (method == "ind") {
    fold <- rep(1:n_folds, each = post_draws)
    df <- data.frame(fold = fold, perfor = fold)
  } else if (method == "avg") {
    fold <- rep(1:n_folds, each = 1)
    df <- data.frame(fold = fold, perfor = fold, lower = fold, upper = fold)
  }
  
  m <- 1
  for (i in 1:n_folds) {
    
    if (method == "ind") {
      for (j in 1:post_draws) {
        real <- kfp$y[pr[[i]]]
        predicted <- kfp$yrep[j, pr[[i]]]
        
        if (metric == "R2") {
          out <- R2(real, predicted)[[1]]
        } else {
          out <- rmse(real, predicted)
        }
        
        df[m, 2] <- out
        m <- m + 1
      }
    } else if (method == "avg") {
      real <- kfp$y[pr[[i]]]
      predicted <- colMeans(kfp$yrep[, pr[[i]]])
      
      if (metric == "R2") {
        out <- R2(real, predicted)
      } else {
        out <- rep(rmse(real, predicted), times = 3)
      }
      
      df[m, 2] <- out[[1]]
      df[m, 3] <- out[[2]]
      df[m, 4] <- out[[3]]
      m <- m + 1
    }
    print(i)
  }
  return(df)
}

recode_func <- function(vector) {
  vector <- if_else(str_detect(vector, "ageE2"), "age^2", vector)
  vector <- if_else(str_detect(vector, "attention_organization"), "attention_diff_organization", as.character(vector))
  vector <- if_else(str_detect(vector, "context_concentrating"), "context_diff_concentrating", as.character(vector))
  vector <- if_else(str_detect(vector, "everything_was_an_effort"), "everything_an_effort", as.character(vector))
  vector <- if_else(str_detect(vector, "nothing_could_cheer"), "sad_no_cheer", as.character(vector))
  vector <- if_else(str_detect(vector, "restless_or_"), "restless_fidgety", as.character(vector))
  vector <- str_remove_all(vector, c("_stress$"))
  vector <- str_remove_all(vector, c("emotions_"))
  vector <- str_remove_all(vector, c("anx_dep_|_Interruptions"))
  
  return(vector)
}