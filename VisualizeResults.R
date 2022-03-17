# -------------------------------
# title: "Visualize Results"
# author: "Zoë Hawks"
# date: "3/10/2022"
# -------------------------------

# -------------------------------
# Part 1: Analysis prep
# -------------------------------
# source functions, read in RDS file, & load libraries
path <- ""
if (path == "") {
  path <- getwd()
} else {
  setwd(path)
}

source("VisualizeResultsFuncs.R")
clean_cor_dat <- readRDS("cleanData.rds")
load_libraries()

if (!dir.exists(paste0(path, "/outputs"))) {
  dir.create(paste0(path, "/outputs"))
} else {
  print("Directory already exists")
}

# -------------------------------
# Part 2: Cog correlations
# -------------------------------
# BP correlations in average performance
for_cor_BP <-
  clean_cor_dat %>% 
  filter(key %in% c("MOT_percent_correct", "GCPT_dprime", "DS_medianRTc", "CRT_medianRTc")) %>%
  unnest(c(data)) %>%
  select(key, value, sitting_id, user_id) %>%
  ungroup() %>%
  mutate(key = str_replace(key, "_", ": ")) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  group_by(user_id) %>%
  summarise_at(vars(contains(":")), mean, na.rm = TRUE) %>%
  ungroup() %>% select(-user_id) %>%
  correlate(use = "pairwise.complete.obs") %>%
  column_to_rownames(., 'term') %>%
  as.matrix(.) %>%
  replace(is.na(.), 1) %>%
  ggcorrplot(., hc.order = TRUE, lab = TRUE, show.diag = FALSE,
             type = "lower",
             colors = c("#6D9EC1", "white", "#E46726")) +
  theme(legend.position = "bottom")

# average WP correlations
# Will throw warning if participant doesn't exhibit variation - cor set to NA, 
# handled in subsequent filter step prior to mean calculation
for_cor_WP <-
  clean_cor_dat %>% 
  filter(key %in% c("MOT_percent_correct", "GCPT_dprime", "DS_medianRTc", "CRT_medianRTc")) %>%
  unnest(c(data)) %>%
  select(key, value, sitting_id, user_id) %>%
  ungroup() %>%
  mutate(key = str_replace(key, "_", ": ")) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  split(.$user_id) %>%
  map(select, -c(sitting_id, user_id)) %>%
  map(cor, use = "pairwise.complete.obs") %>%
  reshape2::melt() %>% select(-L1) %>%
  filter(value != 1 & !is.na(value)) %>%
  ungroup() %>% group_by(Var1, Var2) %>%
  summarise(Mean = mean(value)) %>%
  pivot_wider(names_from = Var2, values_from = Mean) %>%
  column_to_rownames(., 'Var1') %>%
  relocate(`MOT: percent_correct`, .before = `GCPT: dprime`) %>%
  as.matrix(.) %>%
  replace(is.na(.), 1) %>%
  ggcorrplot(., hc.order = TRUE, lab = TRUE, show.diag = FALSE,
             type = "lower",
             colors = c("#6D9EC1", "white", "#E46726")) +
  theme(legend.position = "bottom")

ggsave("outputs/BP_corrmat.tiff", for_cor_BP, width = 8, height = 6, 
       units = "in", bg = "white")
ggsave("outputs/WP_corrmat.tiff", for_cor_WP, width = 8, height = 6, 
       units = "in", bg = "white")

# -------------------------------
# Part 3: Identify predictors in full model
# -------------------------------
# read-in output from BayesianPrediction
files <- list.files("model_results", pattern = "_CVfits.rds$", full.names = TRUE)
rds_list <- map(files, read_rds) %>% bind_rows()

# ICC
ICC_df <-
  rds_list %>%
  select(key, fit_bayes_horse_Int) %>%
  mutate(icc = map(.x = fit_bayes_horse_Int, ~performance::icc(.x))) %>%
  mutate(hoisted = map_dbl(.x = icc, ~.x$ICC_adjusted))

# Model diagnostics
diagnostics <-
  rds_list %>% 
  select(key, fit_bayes_horse_Int, fit_bayes_horse_reTrue) %>%
  gather(key = model, value = value, -key) %>%
  mutate(rhat = map(.x = value, ~(diagnostic_posterior(.x, diagnostic = "Rhat") %>% 
                                    summarise_at(vars(Rhat), .funs = c(min, max))))) %>%
  unnest_wider(c(rhat)) %>%
  rename(min_rhat = fn1, max_rhat = fn2)

# Comparing models based on elpd
# Recommend elpd_diff > 4(se_diff), where larger elpd is better
model_comparisons <-
  rds_list %>% 
  select(key, contains("grouped"), contains("stratefied"), -contains("kfold")) %>%
  gather(key = model, value = value, -key) %>%
  tidyr::extract(model, into = c("model", "method"), regex = "(.*)_([^_]+$)") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(elpd_Int = map_dbl(.x = fit_bayes_horse_IntCV, ~(.x$estimates[1,1])),
         elpd_Full = map_dbl(.x = fit_bayes_horse_reTrueCV, ~(.x$estimates[1,1])),
         se_Int = map_dbl(.x = fit_bayes_horse_IntCV, ~(.x$estimates[1,2])),
         se_Full = map_dbl(.x = fit_bayes_horse_reTrueCV, ~(.x$estimates[1,2])),
         elpd_diff = map2_dbl(.x = fit_bayes_horse_IntCV, fit_bayes_horse_reTrueCV,
                              ~(.x$estimates[1,1] - .y$estimates[1,1])),
         se_diff = map2_dbl(.x = fit_bayes_horse_IntCV, fit_bayes_horse_reTrueCV,
                              ~(.x$estimates[1,2] - .y$estimates[1,2])),
         sig = map2_lgl(.x = elpd_diff, .y = se_diff, ~(abs(.x) > 4*abs(.y))))
model_comparisons %>% 
  select(-method, -contains("fit_bayes")) %>% 
  rename(outcome_var = key) %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = TRUE) %>%
  pack_rows(index = table(model_comparisons$method)) %>%
  cat(., file = "outputs/model_comps.html")
  
# Output table with model results
draws_fixedEffects_unnested <- 
  rds_list %>%
  select(key, fit_bayes_horse_reTrue, fit_bayes_horse_Int) %>% 
  gather(key = model, value = value, -key) %>%
  mutate(model = if_else(str_detect(model, "Int"), "Intercept", "Full"),
         key = map2_chr(.x = key, .y = model, ~paste(.x, .y, sep = "_"))) %>%
  mutate(summarised_tabls = map(.x = value, 
                                ~(broom.mixed::tidy(.x) %>% mutate(term = recode_func(term)) %>%
                                    select(-component)))) %>% 
  select(key, summarised_tabls) %>%
  unnest(c(summarised_tabls)) %>%
  rename(variable = key) %>%
  arrange(variable, estimate)

draws_fixedEffects_unnested %>%
  arrange(variable, estimate) %>%
  select(-variable) %>%
  kbl("html", digits = 3) %>%
  kable_classic(html_font = "Cambria", "basic", "center", full_width = F) %>%
  column_spec(1:7, background = 
                if_else(
                  (draws_fixedEffects_unnested$conf.low < 0 & draws_fixedEffects_unnested$conf.high < 0) |
                    (draws_fixedEffects_unnested$conf.low > 0 & draws_fixedEffects_unnested$conf.high > 0), 
                  "lightgray", "white")) %>%
  pack_rows(index = table(draws_fixedEffects_unnested$variable)) %>%
  cat(., file = "outputs/brms_tables.html")

# Output figure with model results
draws_fixedEffects <- 
  rds_list %>%
  select(key, fit_bayes_horse_reTrue) %>% 
  gather(key = model, value = value, -key) %>%
  mutate(gathered = map(.x = value, ~(gather_draws(.x, `b_.*`, regex = TRUE) )),
         summarised = map(.x = value, ~(gather_draws(.x, `b_.*`, regex = TRUE) %>% 
                                          mean_qi(.width = c(.95, .9, .85, .8)))),
         summarised_abbrev = map(.x = summarised, 
                                 ~filter(.x, (.lower < 0 & .upper < 0) | 
                                           (.lower > 0 & .upper > 0))),
         to_plot = map2(.x = summarised, .y = summarised_abbrev,
                        ~filter(.x, (.variable %in% unique(.y$.variable)))))
plot_fixedEffects <-
  draws_fixedEffects %>%
  mutate(to_plot = map(.x = to_plot, ~(.x %>% filter(.variable != "b_Intercept"))),
         to_plot = map(.x = to_plot, 
                       ~mutate(.x, sig = if_else(((.width >= .95 & .lower < 0 & .upper < 0) |
                                                    (.width >= .95 & .lower > 0 & .upper > 0)), "T", "F"))),
         to_plot = map(.x = to_plot,
                       ~mutate(.x, 
                               .variable = recode_func(.variable),
                               .variable = fct_reorder(.variable, .value))),
         plots = map2(.x = to_plot, .y = key, 
                      ~ggplot(.x, aes(x = .value, y = .variable,
                                      xmin = .lower, xmax = .upper)) +
                        geom_point(shape = 21, size = 5, color = "black",
                                   aes(fill = sig)) + 
                        geom_pointinterval(aes(interval_color = factor(.width)), size = 1) + 
                        labs(title = .y %>% str_replace(., "_", ": ") %>% str_replace(., "_", " "), 
                             x = "Std. beta", y = "", 
                             interval_color = "Interval", fill = "Sig.") +
                        theme_bw() + geom_vline(xintercept = 0, lty = 6) +
                        scale_color_manual(values = c("skyblue", "cornflowerblue", "dodgerblue4", "black"), 
                                           aesthetics = "interval_color") +
                        scale_fill_manual(values = c("transparent", "#66C2A5")) +
                        guides(color = "none", interval_color = guide_legend(override.aes = list(size = 6))) +
                        scale_y_discrete(labels = function(i) str_remove(i, "^[b]_") %>% 
                                           str_remove(., "correct_")) +
                        theme(plot.title = element_text(hjust = 0.5),
                              text = element_text(size = 18))))

# aggregate & arrange to plot
legend <- cowplot::get_legend(
  plot_fixedEffects$plots[[1]] + theme(legend.position = "bottom")
)

comb_full <-
  cowplot::plot_grid(
    cowplot::plot_grid(
      plot_fixedEffects$plots[[1]] + theme(legend.position = "none"), 
      plot_fixedEffects$plots[[2]] + theme(legend.position = "none"), 
      plot_fixedEffects$plots[[3]] + theme(legend.position = "none"), 
      plot_fixedEffects$plots[[4]] + theme(legend.position = "none"), 
      align = "vh", nrow = 2), legend, nrow = 2, rel_heights = c(9, 1))

ggsave("outputs/predictors_full_mod.tiff", comb_full, width = 12, height = 8, 
       units = "in", bg = "white")

# -------------------------------
# Part 4: Predictors across CV runs
# -------------------------------
# evaluating stability of predictors across across folds in grouped and stratified CV
CV_draws <- 
  rds_list %>%
  select(key, fit_bayes_horse_reTrueCV_stratefied, 
         fit_bayes_horse_reTrueCV_grouped, fit_bayes_horse_IntCV_stratefied) %>% 
  gather(key = model, value = value, -key) %>%
  mutate(kfold = map(.x = value, ~unnest_kfold(.x, 5, widths = c(.95, .9, .85, .8))))

CV_draws_new <-
  CV_draws %>%
  mutate(kfold = map(.x = kfold, ~(.x %>% filter(!str_detect(.variable, "Intercept")))),
         summarised_abbrev = map(.x = kfold, 
                                 ~filter(.x, (.lower < 0 & .upper < 0) | 
                                           (.lower > 0 & .upper > 0))),
         to_plot = map2(.x = kfold, .y = summarised_abbrev,
                        ~filter(.x, (.variable %in% unique(.y$.variable)))),
         to_plot = map(.x = to_plot, 
                       ~mutate(.x, sig = if_else(((.width >= .95 & .lower < 0 & .upper < 0) |
                                                    (.width >= .95 & .lower > 0 & .upper > 0)), "T", "F"))),
         to_plot = map(.x = to_plot, 
                       ~mutate(.x, 
                               .variable = recode_func(.variable),
                               .variable = fct_reorder(.variable, .value))),
         to_plot = map(.x = to_plot, 
                       ~mutate(.x, fold = paste0("Fold ", fold))),
         plots = map2(.x = to_plot, .y = key,
                     ~ggplot(.x, aes(x = .value, y = .variable,
                                     xmin = .lower, xmax = .upper)) +
                       geom_point(shape = 21, size = 3, color = "black",
                                  aes(fill = sig)) + 
                       geom_pointinterval(aes(interval_color = factor(.width)), size = 1) + 
                       labs(title = .y %>% str_replace(., "_", ": ") %>% str_replace(., "_", " "), 
                            x = "Std. beta", y = "", 
                            interval_color = "Interval", fill = "Sig.") +
                       theme_bw() + geom_vline(xintercept = 0, lty = 6) +
                       scale_color_manual(values = c("skyblue", "cornflowerblue", "dodgerblue4", "black"), 
                                          aesthetics = "interval_color") +
                       scale_fill_manual(values = c("transparent", "salmon")) +
                       guides(color = "none", interval_color = guide_legend(override.aes = list(size = 6))) +
                       scale_y_discrete(labels = function(i) str_remove(i, "^[b]_") %>% 
                                          str_remove(., "correct_") %>% str_remove(., "_Interruptions")) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~fold, nrow = 1)))

legend2 <- cowplot::get_legend(
  CV_draws_new$plots[[1]] + theme(legend.position = "right")
)

comb2 <- cowplot::plot_grid(
  cowplot::plot_grid(CV_draws_new$plots[[1]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[2]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[3]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[4]] + theme(legend.position = "none"),
                   align = "vh", nrow = 4), 
  legend2, nrow = 1, rel_widths = c(15, 1))

comb3 <- cowplot::plot_grid(
  cowplot::plot_grid(CV_draws_new$plots[[5]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[6]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[7]] + theme(legend.position = "none"), 
                     CV_draws_new$plots[[8]] + theme(legend.position = "none"),
                     align = "vh", nrow = 4), 
  legend2, nrow = 1, rel_widths = c(15, 1))

ggsave("outputs/kfold_stratified.tiff", comb2, width = 15, height = 15, units = "in", bg = "white")
ggsave("outputs/kfold_grouped.tiff", comb3, width = 15, height = 15, units = "in", bg = "white")

# -------------------------------
# Part 5: Overall performance evaluated wrt R2
# -------------------------------
rm(list = setdiff(ls(), c("rds_list", "perf_cal", "recode_func", "R2", "path")))

# one of either "predict" or "fitted"
kfold_method <- "predict"
performance <- 
  rds_list %>%
  select(key, fit_bayes_horse_reTrueCV_stratefied, fit_bayes_horse_reTrueCV_grouped, 
         fit_bayes_horse_IntCV_stratefied, fit_bayes_horse_IntCV_grouped) %>% 
  gather(key = model, value = value, -key) %>%
  mutate(kfp = map(.x = value, ~kfold_predict(.x, method = kfold_method)),
         R2 = map2(.x = kfp, .y = value, ~perf_cal(.y, .x, method = "avg")))

# Compute R2 by model type
performance_summ <- 
  performance %>% 
  mutate(model = str_remove(model, "fit_bayes_horse_")) %>%
  select(-value, -kfp) %>%
  unnest(R2) %>%
  mutate(model = if_else(model == "reTrueCV_stratefied", "(1)", model),
         model = if_else(model == "IntCV_stratefied", "(3)", model),
         model = if_else(model == "reTrueCV_grouped", "(2)", model),
         model = if_else(model == "IntCV_grouped", "(4)", model),
         key = key %>% str_replace(., "_", ": ") %>% str_replace(., "_", " ")) %>%
  group_by(key, model) %>%
  summarise(mean = round(mean(perfor), 4)) %>%
  mutate(mean_chr = paste0("Mean R<sup>2</sup><br>**", mean*100, "%**"))

# Note: in grouped random intercept models, negative correlation messed up sorting LCI vs. UCI 
# When this happened, examined folds and computed LCI, UCI manually
R2_with_CI <- 
  performance %>% 
  mutate(model = str_remove(model, "fit_bayes_horse_")) %>%
  select(-value, -kfp) %>%
  unnest(R2) %>%
  group_by(key, model) %>%
  summarise_at(vars(perfor, lower, upper), .funs = c(mean))

# Visualize R2 by model type, method #1
dodge <- position_dodge(width = 1)  
tmp3 <- 
  performance %>%
  mutate(model = str_remove(model, "fit_bayes_horse_")) %>%
  select(-value, -kfp) %>%
  unnest(R2) %>%
  mutate(model = if_else(model == "reTrueCV_stratefied", "(1)", model),
         model = if_else(model == "IntCV_stratefied", "(3)", model),
         model = if_else(model == "reTrueCV_grouped", "(2)", model),
         model = if_else(model == "IntCV_grouped", "(4)", model),
         key = key %>% str_replace(., "_", ": ") %>% str_replace(., "_", " ")) %>%
  ggplot(., aes(x = fold, y = perfor, fill = factor(fold))) +
  geom_segment(data = performance_summ, 
               aes(y = mean, yend = mean, color = "red4"), x = 0, xend = 6, lty = 2, 
               size = 1, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = .2) +
  geom_point(aes(fill = factor(fold)), position = dodge, shape = 21, size = 2) +
  facet_grid(key~model) +
  geom_richtext(data = performance_summ, 
            aes(label = mean_chr, y = .92, x = .5, hjust = 0), 
            inherit.aes = FALSE, size = 3.5) +
  theme_bw() + labs(x = "Fold", y = expression(R^2),
                    title = "Cross-Validated Prediction", color = "Mean") +
  scale_fill_manual(values = rep("red4", times = 5))  + 
  guides(fill = "none", color = guide_legend(override.aes = list(lty = 1, size = 3))) +
  scale_color_hue(labels = expression(R^2)) +
  theme(text = element_text(size = 12, face = "bold")) +
  ylim(0, 1.02) 

# Visualize R2 by model type, method #2
dodge <- position_dodge(width = .3) 
tmp4 <- 
  performance %>%
  mutate(model = str_remove(model, "fit_bayes_horse_")) %>%
  select(-value, -kfp) %>%
  unnest(R2) %>%
  mutate(model = if_else(model == "reTrueCV_stratefied", "(1)", model),
         model = if_else(model == "IntCV_stratefied", "(3)", model),
         model = if_else(model == "reTrueCV_grouped", "(2)", model),
         model = if_else(model == "IntCV_grouped", "(4)", model),
         key = key %>% str_replace(., "_", ": ") %>% str_replace(., "_", " ")) %>%
  ggplot(., aes(x = model, y = perfor, fill = factor(fold))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = .2) +
  geom_point(aes(fill = factor(fold)), position = dodge, shape = 21, size = 2) +
  geom_point(data = performance_summ,
               aes(y = mean, x = model, color = "red3"), fill = "transparent",
               size = 5, inherit.aes = FALSE, shape = 18) +
  facet_wrap(~key, nrow = 2) +
  geom_richtext(data = performance_summ,
                aes(label = mean_chr, y = rep(c(.9), 16), x = model),
                inherit.aes = FALSE,
                label.padding = unit(c(.15), "lines"),
                size = 4) +
  theme_bw() + labs(x = "", y = "Variance explained",
                    title = "", color = "", fill = "Fold") +
  scale_fill_brewer(palette = "Blues")  + 
  guides(color = guide_legend(override.aes = list(lty = 1, size = 5, order = 2)),
         fill = guide_legend(order = 1)) +
  scale_color_manual(labels = "Mean", values = "red3") +
  theme(text = element_text(size = 18)) +
  ylim(0, 1.05) 

ggsave(paste0("outputs/CV_performance_vis1_", kfold_method, ".tiff"), tmp3, 
       width = 10, height = 8.5, units = "in", bg = "white")
ggsave(paste0("outputs/CV_performance_vis2_", kfold_method, ".tiff"), tmp4, 
       width = 9, height = 7, units = "in", bg = "white")

# ---------------------------------
# Part 6: Visualizing circadian effects
# ---------------------------------
rds_subsetted <- 
  rds_list %>% 
  select(key, data, data_GM, fit_bayes_horse_reTrue) %>%
  mutate(tmp_dat = map(.x = data_GM, ~(.x %>% summarise_if(is.numeric, mean) %>%
                                      mutate(gender = "male", englishPrimary = "yes") %>% 
                                      select(-sin_studytime, -cos_studytime, -age))),
         harm_time = map(.x = data, ~data.frame(cos_time = cos(.x$study_time*2*pi), 
                                                sin_time = sin(.x$study_time*2*pi))),
         time_vec = map(.x = data_GM, ~seq(min(.x$study_time), max(.x$study_time), length.out = 100)),
         summarise_means_sds = map(.x = data, ~(.x %>% select(value, age) %>% 
                                                  summarise_all(list(mean = mean, sd = sd)))))

all_pred <- 
  rds_subsetted %>%
  mutate(expanded = map2(.x = time_vec, .y = harm_time, 
                         ~(expand.grid(time = .x,
                                       age = c(-1, 0, 1))) %>%
                           mutate(cos_studytime = (cos(time*2*pi) - 
                                                     mean(.y$cos_time, na.rm = TRUE))/sd(.y$cos_time, na.rm = TRUE),
                                  sin_studytime = (sin(time*2*pi) - 
                                                     mean(.y$sin_time, na.rm = TRUE))/sd(.y$sin_time, na.rm = TRUE)) %>%
                           select(-time)),
         expanded = map2(.x = expanded, .y = tmp_dat, ~(.x %>% tibble(.y))),
         epred = map2(.x = expanded, .y = fit_bayes_horse_reTrue,
                      ~(.x %>% add_epred_draws(.y, re_formula = NA, scale = "response",
                                               allow_new_levels = TRUE,
                                               ndraws = 1e3))),
         x_vec = map2(.x = epred, .y = harm_time, ~c(asin(.x$sin_studytime*sd(.y$sin_time, na.rm = TRUE) + 
                                                            mean(.y$sin_time, na.rm = TRUE))/(2*pi),
                                                     (pi - asin(.x$sin_studytime*sd(.y$sin_time, na.rm = TRUE) + 
                                                                  mean(.y$sin_time, na.rm = TRUE)))/(2*pi),
                                                     (asin(.x$sin_studytime*sd(.y$sin_time, na.rm = TRUE) + 
                                                             mean(.y$sin_time, na.rm = TRUE)) + 2*pi)/(2*pi),
                                                     (pi - asin(.x$sin_studytime*sd(.y$sin_time, na.rm = TRUE) + 
                                                                  mean(.y$sin_time, na.rm = TRUE)) + 2*pi)/(2*pi))),
         new_dat = map(.x = epred, ~(.x %>% bind_rows(.x) %>% bind_rows(.x) %>% bind_rows(.x))),
         new_dat = map2(.x = new_dat, .y = summarise_means_sds, 
                            ~mutate(.x, 
                                    .epred = round((.epred*.y$value_sd + .y$value_mean), 2),
                                    age = round((age*.y$age_sd + .y$age_mean), 2))), 
         plots = map2(.x = new_dat, .y = x_vec, 
                      ~ggplot(.x, aes(x = .y*24, y = .epred, 4, 
                                      group = interaction(.y, age), color = age %>% as.character())) +
                        geom_smooth(aes(y = .epred, group = paste(age, .draw)), 
                                    alpha = .01, size = .01, se = FALSE, show.legend = FALSE) + 
                        geom_smooth(aes(y = .epred, group = paste(age)), se = FALSE) + 
                        scale_color_brewer(palette = "Set2") +
                        geom_smooth(aes(y = .epred, group = paste(age)), color = "black", lty = "dotted", se = FALSE) + 
                        xlim(9, 22) +
                        theme_bw() ),
         plots = map2(.x = plots, .y = key, 
                      ~(.x + labs(x = "Time (in hours)", y = gsub("^[^_]*_", "", .y), 
                                  color = "Age", title = gsub("_.*","",.y)))))

legend3 <- cowplot::get_legend(
  all_pred$plots[[1]] + theme(legend.position = "bottom")
)

interaction_plots <-
  cowplot::plot_grid(
    cowplot::plot_grid(all_pred$plots[[1]] + theme(legend.position = "none"), 
                       all_pred$plots[[2]] + theme(legend.position = "none"), 
                       all_pred$plots[[3]] + theme(legend.position = "none"), 
                       all_pred$plots[[4]] + theme(legend.position = "none"),
                       align = "vh", nrow = 2), 
    legend3, nrow = 2, rel_heights = c(9, 1))

ggsave("outputs/circadian_plots.tiff", interaction_plots, 
       width = 9, height = 6, units = "in", bg = "white")

