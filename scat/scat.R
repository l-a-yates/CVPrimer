
library(tidyverse); library(tidymodels)
library(pbmcapply); library(randomForest)
library(mgcv); library(RColorBrewer)

rm(list=ls())

MAX_CORES <- 45

#----------
# Data prep
#----------

source("scat_funs.R")
data(scat,package = "caret") # Morphometrics on scat data for bobcat, coyote and gray fox
scat <- prep_data(scat) # N.B. fail = canid (factor-level 1), success = felid (factor-level 2)
run_date <- "2021_05_25"
save_dir <- paste0("fits_",run_date,"/")
if(!dir.exists(save_dir)) dir.create(save_dir)

run <- F # run or load model fits

#scat %>% select(where(is_double)) %>% cor
vars <- names(scat %>% select(-y)); vars




#-------------------------------------------------------------------------------
# Stage 1: use TSS with K-fold CV for variable step selection for logistic model
#-------------------------------------------------------------------------------

set.seed(770) # sample(1e4,1)
cv_data_1 <- vfold_cv(scat,10,50) %>% mutate(thresh_glm = 0.5) # data for 10-fold CV repeated 100 times

if(run){
  varsel_tss_1 <- step_glm(vars, cv_data_1, metric = "tss") # approx 2 mins with 45 cores
  varsel_log_density_1 <- step_glm(vars, cv_data_1, metric = "log_density") 
  saveRDS(varsel_tss_1, paste0(save_dir,"varsel_tss_1_",run_date,".rds"))
  saveRDS(varsel_log_density_1, paste0(save_dir,"varsel_log_density_1_",run_date,".rds"))
} else{
  varsel_tss_1 <- readRDS(paste0(save_dir,"varsel_tss_1_",run_date,".rds"))
  varsel_log_density_1 <- readRDS(paste0(save_dir,"varsel_log_density_1_",run_date,".rds"))
}


tss_data_1 <- varsel_tss_1$metric_step %>% 
  imap(~ .x %>% mutate(model = varsel_tss_1$sel[.y])) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = model, values_from = metric) %>% 
  select(-rep)

ll_data_1 <- varsel_log_density_1$metric_step %>% 
  imap(~ .x %>% mutate(model = varsel_log_density_1$sel[.y])) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = model, values_from = metric) %>% 
  select(-rep)


tss_plot_data_1 <- make_plot_data(tss_data_1, varsel_tss_1$sel)
ll_plot_data_1 <- make_plot_data(ll_data_1, varsel_log_density_1$sel)

tss_plot_data_1 %>% 
  ggplot(aes(model)) +
  geom_point(aes(y = metric), size = 2.5) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod)) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  geom_point(aes(y = metric), shape = 1, size = 5, col = "blue", data = ~ .x %>% filter(model == "mass")) +
  labs(title = "Cross-validation estimates for scat models", x = "Predictor",y = "TSS")

ggsave("plots/scat_modsel_cv_v1.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  



plot_model_comparisons(ll_plot_data_1, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 1: Logistic only, threshold = 0.5, score = log_density")



#------------------------
# 1A Exhasutive selection
#------------------------

options(na.action = "na.fail") # for MuMIn::dredge
0:10 %>% map_dbl(~ choose(10,.x)) %>% sum # 1024 models
forms_all <- glm(y~., data = scat, family = binomial) %>% MuMIn::dredge(evaluate = F) %>% map(as.formula)
tss_dredge <- pbmclapply(forms_all, function(form) fit_confusion_glm(cv_data_1,form)$metric, mc.cores = MAX_CORES) # 3:30 with 45 cores
tss_df <- tss_dredge %>% bind_cols()
tss_summary_all <- make_plot_data(tss_df, levels = names(tss_df)) %>% 
  mutate(dim = forms_all %>% map_dbl(~ .x %>% terms %>% attr("variables") %>% {length(.)-2}))

tss_by_dim <- tss_summary_all %>% 
  group_by(dim) %>% 
  arrange(-metric) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(dim) %>% 
  filter(!dim %in% c(0,9,10))

# colors
display.brewer.pal(3,"Dark2")
p_col <- brewer.pal(3,"Dark2")[3]  

tss_by_dim %>% 
  ggplot(aes(factor(dim, levels = dim))) +
    geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod), col = p_col) +
    geom_point(aes(y = metric), size = 2) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    geom_point(aes(y = metric), shape = 1, size = 6, col = p_col, data = ~ .x %>% filter(dim == 4)) +
    labs(title = "Cross-validation estimates for scat models", x = "Number of parameters",y = "TSS")

#ggsave("plots/scat_modsel_cv_v2A.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  

forms_all[tss_by_dim$model[1:4]]

# top 2% of all models
tss_summary_all %>% 
  arrange(-metric) %>% 
  slice(1:(n()/50)) %>% 
  arrange(dim) %>% 
  group_by(dim) %>% 
  mutate(mlet = c(LETTERS[1:(n())])) %>% 
  mutate(moddim = paste0(dim, mlet)) %>% 
  ggplot(aes(factor(moddim, levels = moddim))) +
  geom_point(aes(y = metric)) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod)) +
  theme_bw() +
  labs(title = "Top 2% of all models", x= "dimension/model", y = "TSS")


tss_summary_all %>% 
  arrange(-metric) %>% 
  slice(1:(n()/10)) %>% 
  arrange(dim) %>% 
  #group_by(dim) %>% 
  ggplot(aes(factor(model, levels = model))) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod), col = p_col) +
  geom_point(aes(y = metric)) +
  geom_point(aes(y = metric), shape = 1, size =6, col = p_col, data = ~ .x %>% filter(metric == max(metric))) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank() ) +
  labs(title = "Top 20% of scat models", x= "Models", y = "TSS") 

#ggsave("plots/scat_modsel_top20_v1.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  



#-------------------------------------------------------------------------------------------------
# Stage 2: Use nested CV to tune hyper-parameters and compare logistic regression to random forest
#-------------------------------------------------------------------------------------------------

# for each outer fold, use inner folds to tune threshold, mtry (rf only), and select variables.
set.seed(7193) # sample(1e4,1)
cv_data_2 <- nested_cv(scat, outside = vfold_cv(v = 10, repeats = 50), inside = vfold_cv(v = 10))

# specify RF hyper-parameters
ntree <- 800
mtry_vec <- 1:(ncol(scat)-1)

# tune models 

top10percent <- tss_summary_all %>% 
  arrange(-metric) %>% 
  slice(1:(n()/10)) %>% 
  pull(model)


if(run){
  rf_tune_values <- tune_rf(cv_data_2, mtry_vec, ntree)  # 1:43 seconds for 50 repeats with 45 cores
  glm_tune_values <- tune_glm(cv_data_2, vars) # 1:38 seconds for 50 repeats with 45 cores
  glm_all_tune_values <- tune_glm_dredge(cv_data_2,forms_all[top10percent])
  saveRDS(glm_tune_values, paste0(save_dir,"glm_tune_values_",run_date,".rds"))
  saveRDS(rf_tune_values, paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  saveRDS(glm_all_tune_values, paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
} else{
  glm_tune_values <- readRDS(paste0(save_dir,"glm_tune_values_",run_date,".rds"))
  rf_tune_values <- readRDS(paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  glm_all_tune_values <- readRDS(paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
}

# add tuned parameters to cv_data
cv_data_2$thresh_rf <- rf_tune_values$threshold
cv_data_2$thresh_glm <- glm_tune_values %>% map_dbl("thresh_glm")
cv_data_2$mtry_rf <- rf_tune_values$mtry
cv_data_2$form_glm <- glm_tune_values %>% map("form_glm")



# fit models
if(run){
  fits_rf_best <- fit_rf(cv_data_2, ntree = ntree, type = "best") # 3 seconds with 45 cores
  fits_rf_all<- fit_rf(cv_data_2, ntree = ntree, type = "all") # 2 seconds with 45 cores
  fits_glm <- fit_confusion_glm(cv_data_2)
  saveRDS(fits_rf_best, paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  saveRDS(fits_rf_all, paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  saveRDS(fits_glm, paste0(save_dir,"fits_glm_",run_date,".rds"))
} else {
  fits_rf_best <- readRDS(paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  fits_rf_all <- readRDS(paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  fits_glm <- readRDS(paste0(save_dir,"fits_glm_",run_date,".rds"))
}



# plot stage 2 results
tss_data_2 <- tibble(glm = fits_glm$metric,
       rf_best = fits_rf_best$metric,
       rf_all = fits_rf_all$metric)

tss_plot_data_2 <- make_plot_data(tss_data_2, names(tss_data_2))

plot_model_comparisons(tss_plot_data_2, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 2: nested CV with tuned hyper-parameter, score = TSS")




# Try a gam
if(F){
v_num <- scat %>% select(where(is.double)) %>% colnames();v_num
v_fac <- setdiff(names(scat)[-1],vars_num); v_fac
varsel_tss_gam_1 <- step_gam(v_num, v_fac, cv_data_1, k = 4)
tss_gam_1 <- varsel_tss_gam_1$metric_step %>% 
  imap(~ .x %>% mutate(model = varsel_log_density_1$sel[.y])) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = model, values_from = metric) %>% 
  select(-rep)

tss_plot_data_gam_1 <- make_plot_data(tss_gam_1, varsel_tss_gam_1$sel)
plot_model_comparisons(tss_plot_data_gam_1, "se_mod") +
  labs(title = "Model comparison - GAM", subtitle = "Stage 1: Logistic only, threshold = 0.5, score = TSS")

#gam_tune_values <- tune_gam(cv_data_2, v_num, v_fac, k = 5)
 
}


