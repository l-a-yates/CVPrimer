
library(tidyverse); library(tidymodels)
library(pbmcapply); library(randomForest)
library(mgcv); library(RColorBrewer)
library(glmnet)
library(dplyr); library(tidyr); library(stringr); library(forcats)

rm(list=ls())
options(na.action = "na.fail") # for MuMIn::dredge

MAX_CORES <- 10

#----------
# Data prep
#----------

source("scat_funs.R")
data(scat,package = "caret") # Morphometrics on scat data for bobcat, coyote and gray fox
scat <- prep_data(scat) # N.B. fail = canid (factor-level 1), success = felid (factor-level 2)
vars <- names(scat %>% select(-y)); vars
scat %>% select(where(is_double)) %>% cor

save_dir <- paste0("files_scat","/")

#---------------------------------------------------------------------------------------
# Part 1: use CM-based metrics with K-fold CV for variable selection for logistic model
#---------------------------------------------------------------------------------------

## 1A: Exhaustive selection (All 1024 models)

run_date <- "2021_09_30" #
#save_dir <- paste0("fits_",run_date,"/")
#if(!dir.exists(save_dir)) dir.create(save_dir)
set.seed(770); cv_data_1 <- vfold_cv(scat,10,50) %>% mutate(thresh = 0.5) # data for 10-fold CV repeated 100 times
run <- F
# specify all 2^10 = 1024 models
forms_all <- glm(y~., data = scat, family = binomial) %>% MuMIn::dredge(evaluate = F) %>% map(as.formula)
dim_all <- forms_all %>% map_dbl(~ .x %>% terms %>% attr("variables") %>% {length(.)-2})

if(run){ # fit all models to all folds and compute metric estimates (220 seconds using 40 cores)
  mcc_dredge <- pbmclapply(forms_all, function(form) fit_confusion_glm(cv_data_1,form, mcc)$metric, mc.cores = MAX_CORES) 
  tss_dredge <- pbmclapply(forms_all, function(form) fit_confusion_glm(cv_data_1,form, tss)$metric, mc.cores = MAX_CORES) 
  ll_dredge <- pbmclapply(forms_all, function(form) fit_ll_glm(cv_data_1,form)$metric, mc.cores = MAX_CORES) 
  saveRDS(mcc_dredge, paste0(save_dir,"mcc_dredge_",run_date,".rds"))
  saveRDS(tss_dredge, paste0(save_dir,"tss_dredge_",run_date,".rds"))
  saveRDS(ll_dredge, paste0(save_dir,"ll_dredge_",run_date,".rds"))
} else {
  mcc_dredge <- readRDS(paste0(save_dir,"mcc_dredge_",run_date,".rds"))
  tss_dredge <- readRDS(paste0(save_dir,"tss_dredge_",run_date,".rds"))
  ll_dredge <- readRDS(paste0(save_dir,"ll_dredge_",run_date,".rds"))
}

metric <- "mcc"
metric_dredge <- get(paste0(metric,"_dredge"))
metric_df <- metric_dredge %>% bind_cols()
metric_summary_all <- make_plot_data(metric_df, levels = names(metric_df)) %>% mutate(dim = dim_all)

# select best scoring model at each dimension
metric_by_dim <- metric_summary_all %>% 
  group_by(dim) %>% 
  arrange(-metric) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(dim) %>% 
  filter(!dim %in% c(0,9,10))

best_dim <- metric_by_dim %>% filter(se_diff == 0) %>% pull(dim)
best_mode_ose_dim <- metric_by_dim %>% filter(metric + se_mod >= max(metric)) %>% filter(dim == min(dim)) %>% pull(dim)
forms_all[metric_by_dim$model[1:5]] # top few models

#display.brewer.pal(3,"Dark2") # colors
p_col <- brewer.pal(3,"Dark2")[3]  

# main plot for scat example in manuscript
metric_by_dim %>% 
  ggplot(aes(factor(dim, levels = dim))) +
    geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod), col = "black") +
    geom_point(aes(y = metric), size = 2) +
    theme_classic() +
    #theme(panel.border = element_blank()) +
    geom_point(aes(y = metric), shape = 1, size = 6, col = "black", data = ~ .x %>% filter(dim == best_mode_ose_dim)) +
    labs(subtitle = "CV estimates for scat models", x = "Number of parameters",y = toupper(metric))

#ggsave("plots/scat_mcc_modsel_2021_10_01.pdf", width = 90, height = 80, units = "mm", device = cairo_pdf()); dev.off()  

# plot for supplementary materials
# top X% of all models
Xpercent = 10

best_mod_ose <- metric_summary_all %>% filter(metric + se_mod >= max(metric, na.rm = T)) %>% 
  filter(dim == min(dim)) %>% filter(metric == max(metric))

metric_summary_all %>% 
  arrange(-metric) %>% 
  slice(1:(n()*(Xpercent/100))) %>% 
  arrange(dim) %>% 
  ggplot(aes(factor(model, levels = model))) +
  geom_linerange(aes(ymin = metric - se_mod, ymax = metric + se_mod, col = factor(dim))) +
  geom_point(aes(y = metric), size =1) +
  geom_point(aes(y = metric), shape = 1, size =5, col = "black", 
             data = ~ .x %>% filter(metric + se_mod >= max(metric)) %>% 
               filter(dim == min(dim))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_color_brewer(name = "Model dimension", palette = "Dark2") +
  theme(panel.border = element_blank(), panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank() ) +
  geom_hline(aes(yintercept = metric), col = "grey20", lty = "dashed", data = ~ .x %>% filter(metric == max(metric))) +
  labs(title = paste0("Top ", Xpercent,"% of scat models"), x= NULL, y = toupper(metric)) 

#ggsave(paste0("plots/scat_modsel_top",Xpercent,"_2021_10_01.pdf"), width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  

# extract and save list of the top 10% of model formulas for stage 2 analysis
top10percent <- metric_summary_all %>%  # use only top 10 percent of models to reduce comp. time (small cheat)
  arrange(-metric) %>% 
  slice(1:(n()/10)) %>% 
  pull(model)

#forms_all[top10percent] %>% saveRDS("fits_2022_01_11/forms_all_top10percent.rds")

#------
## 1B: Penalised regression: lasso, ridge, and elastic net regularisation
#------

run_date <- "2021_09_30" #
save_dir <- paste0("fits_",run_date,"/")
if(!dir.exists(save_dir)) dir.create(save_dir)

set.seed(770) ; cv_data_1 <- vfold_cv(scat,10,50) %>% mutate(thresh = 0.5) # data for 10-fold CV repeated 100 times

# set ranges for regularisation parameter (lambda) values/ranges
log_lambda <- seq(1.5,-4, length.out = 100) # uniformly spaced on the log scale
lambda <- exp(log_lambda)

## tune alpha
## lasso (alpha = 1) is the best performing (MCC) elastic-net model
## ridge performance improves a lot when using alpha = 0.05 instead of 0. 
alpha_vec <-seq(0, 1, 0.025);
if(F){
  metric_alpha <- lapply(alpha_vec, function(a){
    pbmclapply(1:nrow(cv_data_1), function(i) get_cm(cv_data_1$splits[[i]], cv_data_1$id[[i]], cv_data_1$id2[[i]], alpha = a),
               mc.cores = MAX_CORES) %>% bind_rows() %>% as_tibble() %>% 
      group_by(rep, lambda) %>% 
      summarise(cm11 = sum(cm11), cm12 = sum(cm12), cm21 = sum(cm21), cm22 = sum(cm22), .groups = "keep") %>% 
      summarise(metric = (cm11*cm22 - cm12*cm21)/(sqrt((cm11+cm12)*(cm11+cm21)*(cm22+cm12)*(cm22+cm21))), .groups = "drop") %>% 
      group_by(lambda) %>%  
      summarise(metric = mean(metric)) %>% 
      filter(metric == max(metric, na.rm = T)) %>% 
      pull(metric)
  })
  saveRDS(metric_alpha,"fits_scat/metric_alpha.rds")
} else metric_alpha <- readRDS("files_scat/metric_alpha.rds")

tibble(metric = metric_alpha %>% unlist, alpha = alpha_vec) %>% 
  ggplot(aes(alpha,metric)) + geom_line() + theme_classic() +
  labs(subtitle = "Tuning alpha for elastic-net regularisation (ridge = 0, lasso = 1)", y = "MCC")

# set lambda for selected alpha
alpha = 0 # set alpha (1 - lasso, 0 - ridge)
if(alpha == 1) log_lambda <- seq(-1.4,-4, length.out = 100) # lasso
if(alpha <= 0.05) log_lambda <- seq(1.5,-1, length.out = 100)  # (almost) ridge

# compute confusion-matrix entries and aggregate across folds within a given repetition
cm_reg <- pbmclapply(1:nrow(cv_data_1), function(i) get_cm(cv_data_1$splits[[i]], cv_data_1$id[[i]], cv_data_1$id2[[i]], 
                                                           alpha = alpha, lambda = lambda),
                     mc.cores = 10) %>% bind_rows() %>% as_tibble() 

# compute metric for each repetition
metric_reg <- cm_reg %>% 
  group_by(rep, lambda) %>% 
  summarise(metric = matrix(c(sum(cm11),sum(cm21),sum(cm12),sum(cm22)),2,2) %>% mcc)

# find lambda for best score
lambda_best <- metric_reg %>% 
  group_by(lambda) %>%  
  summarise(metric = mean(metric)) %>% 
  filter(metric == max(metric, na.rm = T)) %>% 
  pull(lambda)

# best score
metric_best <- metric_reg %>% 
  filter(lambda == lambda_best) %>%
  arrange(lambda,rep) %>% pull(metric)

# compute summary statistics
metric_plot_data <- metric_reg %>% 
  group_by(lambda) %>% 
  arrange(lambda,rep) %>%
  summarise(se_ose = sd(metric),
            se_diff  = sd(metric - metric_best),
            se_best = sd(metric_best),
            rho_best_m = (se_diff^2 - se_ose^2 - se_best^2)/(-2*se_ose*se_best),
            se_mod = sqrt(1-(rho_best_m))*se_best,
            metric_diff = mean(metric - metric_best),
            metric = mean(metric))


metric_mod_ose <- metric_plot_data %>% 
  filter(metric + se_diff > max(metric, na.rm = T)) %>% 
  filter(lambda == max(lambda))


plot_tuned <- metric_plot_data %>% 
  #filter(metric > 0.3) %>% 
  ggplot(aes(log(lambda))) +
  #geom_linerange(aes(ymax = metric +sd, ymin = metric -sd), col = "grey70", size = 0.8) +
  geom_vline(aes(xintercept = metric_mod_ose$lambda %>% log), col = "grey70", lty = "dashed") +
  geom_vline(aes(xintercept = lambda_best %>% log), col = "grey70", lty = "dashed") +
  geom_linerange(aes(ymax = metric + se_diff, ymin = metric -se_diff), col = "black", size = 0.5, data = metric_mod_ose) +
  geom_line(aes(y = metric), col = "grey30", lty = "solid") + 
  #geom_point(aes(y = metric), size = 1) + 
  geom_point(aes(y = metric), shape = 1, size = 6, data = metric_mod_ose, col = "black") +
  geom_point(aes(y = metric), size = 1, data = metric_mod_ose) +
  #geom_point(aes(y = metric), shape = 1, size = 6, data = ~ .x %>% filter(lambda == lambda_best), col = "black") +
  geom_point(aes(y = metric), size = 1, data = ~ .x %>% filter(lambda == lambda_best), col = "black") +
  labs(subtitle = NULL, x = NULL, y = "MCC") +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank()) 



fit_glmnet <- glmnet(x = scat %>% select(-y) %>% as.data.frame() %>% makeX(),
       y = scat$y,
       family = "binomial",
       lambda = lambda,
       alpha = alpha)


plot_est <- fit_glmnet$beta %>% as.matrix %>% t %>% as_tibble() %>% 
  mutate(lambda = log(lambda), metric = metric_plot_data$metric) %>% 
  pivot_longer(!any_of(c("lambda","metric")), values_to = "estimate", names_to = "coefficient") %>% 
  filter(estimate !=0 | coefficient == "cn") %>% 
  ggplot(aes(lambda,estimate)) + 
  geom_vline(aes(xintercept = metric_mod_ose$lambda %>% log), col = "grey70", lty = "dashed") +
  geom_vline(aes(xintercept = lambda_best %>% log), col = "grey70", lty = "dashed") +
  geom_line(aes(col = coefficient)) +
  #geom_point(aes(col = coefficient), size = 1, data = ~ .x %>% filter(lambda == lambda_best %>% log)) +
  geom_point(aes(col = coefficient), size =1, data = ~ .x %>% filter(lambda == metric_mod_ose$lambda %>% log)) +
  #geom_point(aes(col = coefficient), shape = 1, size =6, data = ~ .x %>% filter(lambda == metric_mod_ose$lambda %>% log)) +
  geom_hline(aes(yintercept = 0), lty = "longdash") +
  labs(y = "Parameter estimates", x = expression(paste("log(",lambda,")"))) +
  theme_classic() +
  theme(legend.position = "none")


# format ridge
if(alpha <= 0.05) {
  x_lim <- c(-1,1.5)
  y_lim <- c(-0.39,NA)
  plot_tuned_ridge <-  plot_tuned + labs(y = NULL, subtitle = "") + xlim(x_lim)
  plot_est_ridge <-  plot_est + labs(y = NULL)+ theme(legend.position = "none") + xlim(x_lim) + ylim(y_lim)
  plot_ridge <- ggpubr::ggarrange(plot_tuned_ridge, plot_est_ridge, ncol = 1, labels = c("B","D")) %>% 
    ggpubr::annotate_figure(top = "Ridge")
}

# format lasso
if(alpha == 1) {
  x_lim <- c(NA,-1.3)
  plot_tuned_lasso <-  plot_tuned + labs(subtitle = "") +
    labs(y = expression("MCC"~phantom(hat(theta)))) +
    xlim(x_lim)
  plot_est_lasso <-  plot_est + theme(legend.position = "none") + 
    labs(y  = expression(widehat(theta))) + 
    xlim(x_lim)
  plot_lasso <- ggpubr::ggarrange(plot_tuned_lasso, plot_est_lasso, ncol = 1, labels = c("A","C")) %>% 
    ggpubr::annotate_figure(top = "LASSO")
}

ggpubr::ggarrange(plot_lasso, plot_ridge, nrow = 1)

#ggsave("plots/scat_reg_tune_2021_10_14.pdf", width = 180, height = 120, 
#      units = "mm", device = cairo_pdf()); dev.off()


#--------------
## Part 2: BPI
#--------------

library(projpred)
library(brms)
library(ggpubr)
library(bayesplot)

run_date <- "2021_12_23" # BPI

# set horse hyperparameters
n <- nrow(scat); n # 91
D <- length(vars); D # 10
p0 <- 3 # prior guess for the number of relevant variables
tau0 <- p0/(D-p0) * 1/sqrt(n) # scale for tau (notice that stan_glm will automatically scale this by sigma)

# fit a reference model
if(F){
  fit.hs <- brm(y ~ ., family=bernoulli(), data=scat,
                  prior=prior(horseshoe(scale_global = tau0, scale_slab = 1), class=b),
                  chains=4, iter=2000, cores = 4,save_pars = save_pars(all = TRUE))
  
  fit.lasso <- brm(y ~ ., family=bernoulli(), data=scat,
                prior= set_prior("lasso(1)"),
                chains=4, iter=2000, cores = 4,save_pars = save_pars(all = TRUE))
  
  fit.ridge <- brm(y ~ ., family=bernoulli(), 
                  prior = prior(normal(0,10)), 
                  data=scat, chains=4, iter=2000, cores = 4,
                  save_pars = save_pars(all = TRUE))
  
  fit.flat <- brm(y ~ ., family=bernoulli(), 
                   data=scat, chains=4, iter=2000, cores = 4,
                   save_pars = save_pars(all = TRUE))
  
  saveRDS(fit.hs, paste0(save_dir,"fit.hs.rds"))
  saveRDS(fit.lasso, paste0(save_dir,"fit.lasso.rds"))
  saveRDS(fit.ridge, paste0(save_dir,"fit.ridge.rds"))
  saveRDS(fit.flat, paste0(save_dir,"fit.flat.rds"))
} else {
  fit.hs <- readRDS(paste0(save_dir,"fit.hs.rds"))
  fit.lasso <- readRDS(paste0(save_dir,"fit.lasso.rds"))
  fit.ridge <- readRDS(paste0(save_dir,"fit.ridge.rds"))
  fit.flat <- readRDS(paste0(save_dir,"fit.flat.rds"))
}


# loo check
fits <- list(hs = fit.hs, lasso = fit.lasso, ridge = fit.ridge, flat = fit.flat)

if(F){
  future::plan("multisession", workers = 10)
  fits.loo <- fits %>% map(loo, reloo = T, future = T) # hs (best) -> lasso -> weak -> flat (worst)
  saveRDS(fits.loo, paste0(save_dir,"fits.loo.rds"))
  
} else {
  fits.loo <- readRDS(paste0(save_dir,"fits.loo.rds"))
}

fits.loo %>% loo_compare %>% print(simplify = F)

# plot posteriors of hs model
plot_post <- function(fit){
  fit %>% as_tibble %>% select(-starts_with("l")) %>% mcmc_areas(area_method = "scaled", prob_outer = 0.98) +
    xlim(c(-2.8,2)) +
    theme_classic()
}

ggarrange(
  fit.ridge %>% plot_post + labs(subtitle = "Weakly informative"),
  fit.lasso %>% plot_post + labs(subtitle = "Lasso"),
  fit.hs %>% plot_post + labs(subtitle = "Horseshoe"),
  ncol = 1,
  labels = "AUTO"
)
# ggsave("plots/scat_reg_post.pdf", plot = ., width = 150, height = 200, units = "mm")


# varsel and projection
if(F){
  vs <- cv_varsel(fit.hs, cv_method = "LOO", method = "forward")
  saveRDS(vs,paste0(save_dir,"vs.rds"))
} else vs <- readRDS(paste0(save_dir,"vs.rds"))

# plot varsel
vs_plot <- 
  vs %>% plot(stats = c("elpd"), deltas = T) + 
  theme_classic() +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  labs(y = expression(Delta*"ELPD"))


vs_plot
solution_terms(vs)

# project onto 1 parameter submodel
proj <- project(vs, nterms = 1, ndraws = 4000)

# fit single model -- no projection
if(F){
  fit.cn <- brm(y ~ cn, family=bernoulli(), data=scat,chains=4, iter=2000, cores = 4)
  saveRDS(fit.cn, paste0(save_dir,"fit.cn.rds"))
} else fit.cn <- readRDS(paste0(save_dir,"fit.cn.rds"))

plot.ref <- as.data.frame(fit.hs) %>% select(Intercept = b_Intercept, cn = b_cn) %>% 
  mcmc_areas() + labs(subtitle = "Reference") + xlim(c(-3.5,2))

plot.cn <- as.data.frame(fit.cn) %>% select(Intercept = b_Intercept, cn = b_cn) %>% 
  mcmc_areas() + labs(subtitle = "CN only") + xlim(c(-3.5,2))

plot.proj <- mcmc_areas(as.matrix(proj))+ labs(subtitle = "Projected") + xlim(c(-3.5,2))

ggarrange(plot.ref, plot.proj, plot.cn, nrow = 3, ncol= 1)

# compare posteriors
cn_post <- 
  tibble(ref = fit.hs %>% as_tibble() %>% pull(b_cn),
       proj = proj %>% as.matrix %>% {.[,"cn"]},
       `non-proj` = fit.cn %>% as_tibble() %>% pull(b_cn)) %>% 
  mcmc_areas() +
  labs(x = expression(theta["carbon-nitrogen ratio"])) +
  theme_classic()


ggarrange(vs_plot + labs(subtitle = "Projective step selection"), 
          cn_post + labs(subtitle = "Posterior density", y = "Model"), 
          labels = "AUTO",
          ncol = 2)

#ggsave("plots/scat_BPI.pdf", width = 7.3, height = 3)



#-------------------------------------------------------------------------------------------------
# Part 3: nested CV to compare logistic regression to random forest
#-------------------------------------------------------------------------------------------------

run_date <- "2022_01_11" # nested analysis with MCC

# for each outer fold, use inner folds to tune threshold, mtry (rf only), and select variables.
set.seed(7193) # sample(1e4,1)
cv_data_2 <- nested_cv(scat, outside = vfold_cv(v = 10, repeats = 50), inside = vfold_cv(v = 10))
vars <- names(scat %>% select(-y)); vars
forms_top10 <- readRDS(paste0(save_dir,"forms_all_top10percent.rds"))

# specify RF hyper-parameters
ntree <- 500
mtry_vec <- 1:(ncol(scat)-1)

# tune models 
if(run){
  rf_tune_values <- tune_rf(cv_data_2, mtry_vec, ntree, metric = mcc)  # 1:17 seconds for 50 repeats with 40 cores
  glm_step_tune_values <- tune_glm_step(cv_data_2, vars, metric = mcc) # 1:40 seconds for 50 repeats with 40 cores
  glm_all_tune_values <- tune_glm_all(cv_data_2,forms_top10, metric = mcc) # 3:20 seconds for top 10% of model with 40 cores
  saveRDS(glm_step_tune_values, paste0(save_dir,"glm_step_tune_values_",run_date,".rds"))
  saveRDS(rf_tune_values, paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  saveRDS(glm_all_tune_values, paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
} else{
  glm_step_tune_values <- readRDS(paste0(save_dir,"glm_step_tune_values_",run_date,".rds"))
  rf_tune_values <- readRDS(paste0(save_dir,"rf_tune_values_",run_date,".rds"))
  glm_all_tune_values <- readRDS(paste0(save_dir,"glm_all_tune_values_",run_date,".rds"))
}


# add tuned parameters to cv_data
cv_data_2$thresh_rf <- rf_tune_values$threshold
cv_data_2$thresh_glm_step <- glm_step_tune_values$threshold
cv_data_2$thresh_glm_all <- glm_all_tune_values$threshold
cv_data_2$mtry_rf <- rf_tune_values$mtry
cv_data_2$form_glm_step <- glm_step_tune_values$form
cv_data_2$form_glm_all <- glm_all_tune_values$form

cv_data_2 %>% pull(mtry_rf) %>% table
cv_data_2 %>% pull(form_glm_step) %>% str_split(" +") %>% map_dbl(~ .x %in% vars %>% sum) %>% table
#readRDS("files_scat/cv_data_2.rds")

# fit models
if(run){
  fits_rf_best <- fit_rf(cv_data_2, ntree = ntree, type = "best", metric = mcc) # 3 seconds with 40 cores
  fits_rf_all<- fit_rf(cv_data_2, ntree = ntree, type = "all", metric = mcc) # 2 seconds with 40 cores
  fits_glm_step <- fit_confusion_glm(cv_data_2 %>% mutate(form_glm = form_glm_step, thresh = thresh_glm_step), metric = mcc)
  fits_glm_all <- fit_confusion_glm(cv_data_2 %>% mutate(form_glm = form_glm_all, thresh = thresh_glm_all), metric = mcc)
  saveRDS(fits_rf_best, paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  saveRDS(fits_rf_all, paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  saveRDS(fits_glm_step, paste0(save_dir,"fits_glm_step_",run_date,".rds"))
  saveRDS(fits_glm_all, paste0(save_dir,"fits_glm_all_",run_date,".rds"))
} else {
  fits_rf_best <- readRDS(paste0(save_dir,"fits_rf_best_",run_date,".rds"))
  fits_rf_all <- readRDS(paste0(save_dir,"fits_rf_all_",run_date,".rds"))
  fits_glm_step <- readRDS(paste0(save_dir,"fits_glm_step_",run_date,".rds"))
  fits_glm_all <- readRDS(paste0(save_dir,"fits_glm_all_",run_date,".rds"))
}


# plot stage 2 results
tss_data_2 <- tibble(glm_step = fits_glm_step$metric,
                     glm_all = fits_glm_all$metric,
                     rf_best = fits_rf_best$metric,
                     rf_all = fits_rf_all$metric)

tss_plot_data_2 <- make_plot_data(tss_data_2, names(tss_data_2))

#  model    metric     se se_diff se_mod
#  <fct>     <dbl>  <dbl>   <dbl>  <dbl>
#1 glm_step  0.250 0.0677  0.0969 0.0639
#2 glm_all   0.250 0.0645  0.0946 0.0641
#3 rf_best   0.375 0.0594  0      0     
#4 rf_all    0.355 0.0639  0.0603 0.0410

plot_model_comparisons(tss_plot_data_2, "se_mod") +
  labs(title = "Model comparison", subtitle = "Stage 2: nested CV with tuned hyper-parameter, score = TSS")



#-------------------
# bias-variance plot
#-------------------

arr1<- arrow(ends = "both", type = "closed", length = unit(2, "mm"))
arr2<- arrow(type = "open", length = unit(2, "mm"))
tibble(x = seq(0,5,0.01),
       bias = exp(-x) + 0.2,
       var = exp(1.2*x)/exp(1.2*5)+0.2,
       tot = bias + var + 0.08) %>% 
  ggplot(aes(x)) +
  geom_line(aes(y = bias), col = brewer.pal(8,"Dark2")[1], arrow = arr1) +
  geom_line(aes(y = var), col = brewer.pal(8,"Dark2")[2] , arrow = arr1) +
  geom_line(aes(y = tot), arrow = arr1) +
  annotate(geom = "text", col = brewer.pal(8,"Dark2")[1], x = 0.13, y = 0.8, label = expression(bias^2)) +
  annotate(geom = "text", col = brewer.pal(8,"Dark2")[2], x = 0.3, y = 0.3, label = "variance") +
  annotate(geom = "text", col = "black", x = 1.1, y = 1.1, label = "total error") +
  ylim(c(0.13,1.2)) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(arrow = arr2)) +
  labs(x = "Complexity", y = "Error") 

#ggsave("plots/bias_variance_plot.pdf", width = 4, height  = 3)

