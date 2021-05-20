rm(list=ls())

library(tidyverse); library(tidymodels)
library(pbmcapply); library(randomForest)


#----------
# Data prep
#----------

source("scat_funs.R")
data(scat,package = "caret") # Morphometrics on scat data for bobcat, coyote and gray fox
scat <- prep_data(scat) # N.B. Success = canid, Failure = felid
scat %>% group_by(number) %>% summarise(canid = sum(y=="canid"), felid = sum(y=="felid"))
scat
scat <- scat %>% mutate(across("number",as.numeric))

#scat <- scat %>% mutate(number = pmin(number,5)) 
#scat <- scat %>% mutate(number = pmin(number,5) %>% factor(levels = c("1","2","3","4","5"), labels = c("1","2","3","4","5+")))
#scat <- scat %>% select(-number)

set.seed(3218) # sample(1e4,1)
cv_data <- vfold_cv(scat,10,100,strata = number, bins = 5)

vars_num <- scat %>% select(where(is.double)) %>% colnames();vars_num
vars_fac <- setdiff(names(scat)[-1],vars_num) 
make_form <- function(v_num,v_fac, k = 7) paste("y ~", paste("s(",v_num,", k = ", k , ")", collapse = " + "), " + ", paste(v_fac, collapse = " +")) %>% as.formula()
make_form(vars_num, vars_fac)
mgcv::gam(make_form(vars_num, vars_fac), family = binomial(), data = scat) %>% plot
run_date <- "2021_05_18"
save_dir <- paste0("fits_",run_date,"/")


#----------------------------------------
# tune mtry and threshold using TSS. 
# Apply to all CV splits & add to cv_data
#-----------------------------------------

mtry_vec <- 1:(ncol(scat)-1)
tune_mtry <- cv_data$splits %>% map(analysis) %>% pbmclapply(tune_tss_mtry, mtry_vec = mtry_vec, mc.cores = 40)
#saveRDS(tune_mtry,paste0(save_dir,"tune_mtry_2021_05_18.rds"))
#tune_mtry <- readRDS(paste0(save_dir,"tune_mtry_",run_date,".rds"))
thresh_mtry <- tune_mtry %>% map_dfr(~ .x %>% filter(tss == min(tss)) %>% slice(1) %>% as.vector)
cv_data$thresh <- thresh_mtry$threshold
cv_data$mtry <- thresh_mtry$mtry


tss_rf <- fit_confusion_rf(cv_data)
#saveRDS(tss_rf,paste0(save_dir,"tss_rf_2021_05_18.rds"))
#tss_rf <- readRDS(paste0(save_dir,"tss_rf_2021_05_18.rds"))

#-----------------------------------
# step selection for logistic models
#-----------------------------------

if(T){
  
  vars <- names(scat %>% select(-y)); vars
  sel <- c()
  tss_step <- list()
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
  
  for(step in 1:length(vars)){
    message(paste("Step:",step))
    newfits = pbmclapply(setdiff(vars,sel), function(var){
      fit_confusion_glm(make_formula(c(sel,var)),cv_data)
    }, mc.cores = length(setdiff(vars,sel)))
    best = newfits %>% map("tss") %>% map_dbl(mean) %>% which.max()
    sel[step] <- setdiff(vars,sel)[best]
    message(paste("Selected:",sel[step]))
    tss_step[[step]] <- newfits[[best]]
  }
  varsel <- list(sel = sel, tss_step = tss_step)
  #saveRDS(varsel,paste0(save_dir,"varsel_2021_05_18.rds"))
}

#varsel <- readRDS(paste0(save_dir,"varsel_2021_05_18.rds"))

#------------------
# Plot tss results
#------------------

tss_data <- varsel$tss_step %>% 
  imap(~ .x %>% mutate(model = varsel$sel[.y])) %>% 
  bind_rows() %>% 
  bind_rows(tss_rf %>% mutate(model = "rf")) %>% 
  mutate(model = factor(model, levels = c(varsel$sel,"rf")))

tss_data %>% 
  group_by(model) %>% 
  summarise(TSS = mean(tss), se = sd(tss)/sqrt(n())) %>% 
  ggplot(aes(model)) +
  geom_point(aes(y = TSS)) +
  geom_linerange(aes(ymin = TSS - se, ymax = TSS + se)) +
  theme_bw() +
  labs(title = "with number")
  


