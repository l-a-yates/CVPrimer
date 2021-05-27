
library(tidyverse); library(tidymodels)
library(pbmcapply); library(randomForest)
library(mgcv); library(RColorBrewer)

rm(list=ls())
options(na.action = "na.fail") # for MuMIn::dredge

MAX_CORES <- 20

#----------
# Data prep
#----------

source("scat_funs.R")
data(scat,package = "caret") # Morphometrics on scat data for bobcat, coyote and gray fox
scat <- prep_data(scat) # N.B. fail = canid (factor-level 1), success = felid (factor-level 2)
run_date <- "2021_05_26"
save_dir <- paste0("fits_",run_date,"/")
if(!dir.exists(save_dir)) dir.create(save_dir)

run <- T # run or load model fits

#scat %>% select(where(is_double)) %>% cor
vars <- names(scat %>% select(-y)); vars


#-------------------------------------------------------------------------------
# Stage 1: use TSS with K-fold CV for variable selection for logistic model
#-------------------------------------------------------------------------------

set.seed(770) # sample(1e4,1)
cv_data_1 <- vfold_cv(scat,10,50) %>% mutate(thresh = 0.5) # data for 10-fold CV repeated 100 times

cv_data_1
# Stage 1A : Step selection (TSS evaluation at each step)

if(run){
  varsel_tss_1 <- step_glm(vars, cv_data_1, metric = "tss") # approx 80 seconds with 45 cores
  varsel_log_density_1 <- step_glm(vars, cv_data_1, metric = "log_density") # approx 50 seconds with 45 cores
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

# main plot for step selection
tss_plot_data_1 %>% 
  ggplot(aes(model)) +
  geom_point(aes(y = metric), size = 2.5) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod)) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  geom_point(aes(y = metric), shape = 1, size = 5, col = "blue", data = ~ .x %>% filter(model == "mass")) +
  labs(title = "Cross-validation estimates for scat models", x = "Predictor",y = "TSS")

#ggsave("plots/scat_modsel_cv_v1.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  

# plot for supplementary materials
plot_model_comparisons(ll_plot_data_1, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 1: Logistic only, threshold = 0.5, score = log_density")



# Stage 1B: Exhaustive selection (All 1024 models)

# specify all 2^10 = 1024 models
forms_all <- glm(y~., data = scat, family = binomial) %>% MuMIn::dredge(evaluate = F) %>% map(as.formula)
dim_all <- forms_all %>% map_dbl(~ .x %>% terms %>% attr("variables") %>% {length(.)-2})

if(run){ # fit all models to all folds and compute TSS estimates (400 seconds using 45 cores)
  tss_dredge <- pbmclapply(forms_all, function(form) fit_confusion_glm(cv_data_1,form)$metric, mc.cores = MAX_CORES) 
  saveRDS(tss_dredge, paste0(save_dir,"tss_dredge_",run_date,".rds"))
} else tss_dredge <- readRDS(paste0(save_dir,"tss_dredge_",run_date,".rds"))


tss_df <- tss_dredge %>% bind_cols()
tss_summary_all <- make_plot_data(tss_df, levels = names(tss_df)) %>% 
  mutate(dim = dim_all)

# select best scoring model at each dimension
tss_by_dim <- tss_summary_all %>% 
  group_by(dim) %>% 
  arrange(-metric) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(dim) %>% 
  filter(!dim %in% c(0,9,10))

forms_all[tss_by_dim$model[1:5]] # top few models

#display.brewer.pal(3,"Dark2") # colors
p_col <- brewer.pal(3,"Dark2")[3]  

# main plot for scat example in manuscript
tss_by_dim %>% 
  ggplot(aes(factor(dim, levels = dim))) +
    geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod), col = p_col) +
    geom_point(aes(y = metric), size = 2) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    geom_point(aes(y = metric), shape = 1, size = 6, col = p_col, data = ~ .x %>% filter(dim == 4)) +
    labs(title = "Cross-validation estimates for scat models", x = "Number of parameters",y = "TSS")

#ggsave("plots/scat_modsel_cv_v2A.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  

# plot for supplementary materials
# top X% of all models
Xpercent = 5
tss_summary_all %>% 
  arrange(-metric) %>% 
  slice(1:(n()*(Xpercent/100))) %>% 
  arrange(dim) %>% 
  ggplot(aes(factor(model, levels = model))) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod), col = p_col) +
  geom_point(aes(y = metric)) +
  geom_point(aes(y = metric), shape = 1, size =6, col = p_col, data = ~ .x %>% filter(metric == max(metric))) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank() ) +
  labs(title = paste0("Top ", Xpercent,"% of scat models"), x= "Models", y = "TSS") 

#ggsave(paste0("plots/scat_modsel_top",Xpercent,"_v1.pdf"), width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  

# extract list of top 10% of model names for stage 2 analysis
top10percent <- tss_summary_all %>%  # use only top 10 percent of models to reduce comp. time (small cheat)
  arrange(-metric) %>% 
  slice(1:(n()/10)) %>% 
  pull(model)



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
if(run){
  rf_tune_values <- tune_rf(cv_data_2, mtry_vec, ntree)  # 1:43 seconds for 50 repeats with 45 cores
  glm_step_tune_values <- tune_glm_step(cv_data_2, vars) # 1:38 seconds for 50 repeats with 45 cores
  glm_all_tune_values <- tune_glm_all(cv_data_2,forms_all[top10percent])
  saveRDS(glm_tune_values, paste0(save_dir,"glm_step_tune_values_",run_date,".rds"))
  saveRDS(rf_tune_values, paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  saveRDS(glm_all_tune_values, paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
} else{
  glm_step_tune_values <- readRDS(paste0(save_dir,"glm_step_tune_values_",run_date,".rds"))
  rf_tune_values <- readRDS(paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  glm_all_tune_values <- readRDS(paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
}
rf_tune_values
glm_step_tune_values
# add tuned parameters to cv_data
cv_data_2$thresh_rf <- rf_tune_values$threshold
cv_data_2$thresh_glm_step <- glm_step_tune_values$threshold
cv_data_2$thresh_glm_all <- glm_all_tune_values$threshold
cv_data_2$mtry_rf <- rf_tune_values$mtry
cv_data_2$form_glm_step <- glm_step_tune_values$form
cv_data_2$form_glm_all <- glm_all_tune_values$form


# fit models
if(run){
  fits_rf_best <- fit_rf(cv_data_2, ntree = ntree, type = "best") # 3 seconds with 45 cores
  fits_rf_all<- fit_rf(cv_data_2, ntree = ntree, type = "all") # 2 seconds with 45 cores
  fits_glm_step <- fit_confusion_glm(cv_data_2 %>% mutate(form_glm = form_glm_step, thresh = thresh_glm_step))
  fits_glm_all <- fit_confusion_glm(cv_data_2 %>% mutate(form_glm = form_glm_all, thresh = thresh_glm_all))
  saveRDS(fits_rf_best, paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  saveRDS(fits_rf_all, paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  saveRDS(fits_glm_step, paste0(save_dir,"fits_glm_step_",run_date,".rds"))
  saveRDS(fits_glm_all, paste0(save_dir,"fits_glm_all_",run_date,".rds"))
} else {
  fits_rf_best <- readRDS(paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  fits_rf_all <- readRDS(paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  fits_glm_step <- readRDS(paste0(save_dir,"fits_glm_step",run_date,".rds"))
  fits_glm_all <- readRDS(paste0(save_dir,"fits_glm_",run_date,".rds"))
}


# plot stage 2 results
tss_data_2 <- tibble(glm_step = fits_glm_step$metric,
                     glm_all = fits_glm_all$metric,
                     rf_best = fits_rf_best$metric,
                     rf_all = fits_rf_all$metric)


tss_plot_data_2 <- make_plot_data(tss_data_2, names(tss_data_2))

plot_model_comparisons(tss_plot_data_2, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 2: nested CV with tuned hyper-parameter, score = TSS")

# all-version
tss_data_2_all <- tibble(glm_all = fits_glm_all$metric,
                     rf_all = fits_rf_all$metric)


tss_plot_data_2_all <- make_plot_data(tss_data_2_all, names(tss_data_2_all))

plot_model_comparisons(tss_plot_data_2_all, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 2: nested CV with tuned hyper-parameter",
       y = "TSS")

tss_plot_data_2_all



