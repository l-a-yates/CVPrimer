

library(tidyverse); library(tidymodels)
library(pbmcapply); library(randomForest)
library(mgcv)

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

#scat %>% select(where(is_double)) %>% cor
vars <- names(scat %>% select(-y)); vars

#-------------------------------------------------------------------------------
# Stage 1: use TSS with K-fold CV for variable step selection for logistic model
#-------------------------------------------------------------------------------

set.seed(770) # sample(1e4,1)
cv_data_1 <- vfold_cv(scat,10,50) %>% mutate(thresh_glm = 0.5) # data for 10-fold CV repeated 100 times
varsel_tss_1 <- step_glm(vars, cv_data_1, metric = "tss") # approx 2 mins with 45 cores
varsel_log_density_1 <- step_glm(vars, cv_data_1, metric = "log_density") 


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

plot_model_comparisons(tss_plot_data_1, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 1: Logistic only, threshold = 0.5, score = TSS")
plot_model_comparisons(ll_plot_data_1, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 1: Logistic only, threshold = 0.5, score = log_density")



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
rf_tune_values <- tune_rf(cv_data_2, mtry_vec, ntree)  # 1:43 seconds for 50 repeats with 45 cores
glm_tune_values <- tune_glm(cv_data_2) # 1:38 seconds for 50 repeats with 45 cores

# add tuned parameters to cv_data
cv_data_2$thresh_rf <- rf_tune_values$threshold
cv_data_2$thresh_glm <- glm_tune_values %>% map_dbl("thresh_glm")
cv_data_2$mtry_rf <- rf_tune_values$mtry
cv_data_2$form_glm <- glm_tune_values %>% map("form_glm")

# fit models
fits_rf_best <- fit_rf(cv_data_2, ntree = ntree, type = "best") # 3 seconds with 45 cores
fits_rf_all<- fit_rf(cv_data_2, ntree = ntree, type = "all") # 2 seconds with 45 cores
fits_glm <- fit_confusion_glm(cv_data_2)

# plot stage 2 results
tss_data_2 <- tibble(glm = fits_glm$metric,
       rf_best = fits_rf_best$metric,
       rf_all = fits_rf_all$metric)

tss_plot_data_2 <- make_plot_data(tss_data_2, names(tss_data_2))

plot_model_comparisons(tss_plot_data_2, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 2: nested CV with tuned hyper-parameter, score = TSS")


## a quick gam
if(F){
  vars_num <- scat %>% select(where(is.double)) %>% colnames();vars_num
  vars_fac <- setdiff(names(scat)[-1],vars_num) 
  make_form <- function(v_num,v_fac, k = 7) paste("y ~", paste("s(",v_num,", k = ", k , ")", collapse = " + "), " + ", paste(v_fac, collapse = " +")) %>% as.formula()
  make_form(vars_num, vars_fac)
  fit <- mgcv::gam(make_form(vars_num, vars_fac), family = binomial(), data = scat);fit
  fit %>% AIC
}

