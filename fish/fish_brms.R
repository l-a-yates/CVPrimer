#----------------------------------------------------------------------------
#
# Cross validation for model selection: a primer with examples from ecology.
#
# Code for example 2: Pinfish growth models
#
# Authors: L.A.Yates, ...
#
#----------------------------------------------------------------------------

library(tidyverse)
library(future)
library(fishmethods)
library(FlexParamCurve)
library(RColorBrewer)
library(brms)

rm(list=ls())

#source("fish_functions_nlme.R")
options(brms.file_refit = "on_change", auto_write = TRUE)

run_date <- "2021_10_06"
#sv_dir <- paste0("fits_",run_date)
#if(!dir.exists(sv_dir)) dir.create(sv_dir)

sv_dir <- paste0("/media/layates/CVprimer/fits_",run_date)
if(!dir.exists(sv_dir)) dir.create(sv_dir) else message(paste0("Saving results to ",sv_dir,"/"))

#-------------------------
# Load and prepare data
#-------------------------

data(pinfish); pinfish
names(pinfish)[1] <- "haul"

# remove samples with undetermined sex
fishData <- pinfish[pinfish$sex != 0,]
fishData$sex <- factor(fishData$sex)

# restrict minimum haul size
fishData %>% group_by(haul) %>% summarise(n = n()) %>% pull(n) %>% table() # examine haul sizes
min.size <- 5
fishData <- fishData %>% group_by(haul) %>% filter(n()>= min.size) %>% ungroup
fishData <- fishData %>% filter(haul != "TBD930066") # remove outlier
fishData$haul <-  factor(fishData$haul)

fishData
# MLE fits for initial values and priors
#f <- fit.re(fishData)
#f %>% map_dbl(AIC) %>% {.-min(.)}



  
  nlf.vB <- sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0))))
  nlf.G <- sl~ (152 + Linf) * exp(-exp(-(1.28 + K) * (age - t0)))
  nlf.log <- sl~ (150 + Linf) * exp(-exp(-(1.28 + K) * (age - t0)))

  make_form <- function(fun, K, t, L) {
    if(K) {f.K <- K ~ sex} else {f.K <- K ~ 1}
    if(t) {f.t <- t0 ~ sex} else {f.t <- t0 ~ 1}
    if(L) {f.L <- Linf ~ sex + (1|haul)} else {f.L <- Linf ~ 1 + (1|haul)}
    bf(get(paste0("nlf.",fun)), f.K, f.t, f.L, nl = T)
  }
  
  make_prior <- function(K, t, L){
    p <-  
      prior(normal(0, 1), nlpar = "K") +  
      prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
      prior(normal(0, 1), nlpar = "t0") +
      prior(normal(0, 10), nlpar = "Linf")
    if(K) p <-  p + prior(normal(0, 0.3), nlpar = "K", coef = "sex2") 
    if(t) p <-  p + prior(normal(0, 0.3), nlpar = "t0", coef = "sex2") 
    if(L) p <-  p + prior(normal(0, 0.3), nlpar = "Linf", coef = "sex2") 
    return(p)
  }
  
  make_file_name <- function(x){
    s_type <- paste0(c("K","L","t")[x[1,2:4]==1], collapse = "")
    if(s_type=="") s_type <- "0"
    paste0(sv_dir,paste(c("/m",x[1],s_type), collapse = "."))
  }
  


  make_form("log",T,T,T)
  make_prior(T,T,F)
  
  models_grid <- expand.grid(fun = c("vB","G","log"), K = c(0,1), L = c(0,1), t  = c(0,1), stringsAsFactors = F)
  
  
  asasas
  
  for(i in 1:nrow(models_grid)){
    print(paste("Fitting model", i,make_file_name(models_grid[i,]))) 
    #brm(formula = models_grid[i,] %>% with(make_form(fun,K,t,L)),
    #    prior = models_grid[i,] %>% with(make_prior(fun,K,t,L)),
    #    data = fishData,
    #    cores = 4, 
    #    iter = 4000, 
    #    file = make_file_name(models_grid[i,]))
  }

  
  
  #----------
  # vB models
  #----------
  if(F){
    
  f.vB.0 <-  bf(sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
  f.vB.K <-  bf(sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ sex, t0 ~ 1, Linf ~ 1 + 1|haul, nl = TRUE)
  f.vB.t <-  bf(sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ sex, Linf ~ 1 + 1|haul, nl = TRUE)
  f.vB.Kt <-  bf(sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ sex, t0 ~ sex, Linf ~ 1 + 1|haul, nl = TRUE)
  f.vB.L <-  bf(sl ~ (152 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ 1, Linf ~ sex + (1|haul), nl = TRUE)
  # N.B. vB.L does not improve over vB.0 or vB.t
  
  pr.vB <- 
    prior(normal(0, 1), nlpar = "K") +  
    prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
    prior(normal(0, 1), nlpar = "t0") +
    prior(normal(0, 10), nlpar = "Linf")
  
  pr.vB.K <- pr.vB + prior(normal(0, 0.3), nlpar = "K", coef = "sex2") 
  pr.vB.t <- pr.vB + prior(normal(0, 0.3), nlpar = "t0", coef = "sex2") 
  pr.vB.Kt <- pr.vB.K + prior(normal(0, 0.3), nlpar = "t0", coef = "sex2") 
  pr.vB.L <- pr.vB + prior(normal(0, 0.3), nlpar = "Linf", coef = "sex2") 
  
  m.vB.0 <- brm(f.vB.0, data = fishData, prior = pr.vB, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.0"))
  m.vB.K <- brm(f.vB.K, data = fishData, prior = pr.vB.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.K"))
  m.vB.t <- brm(f.vB.t, data = fishData, prior = pr.vB.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.t"))
  m.vB.Kt <- brm(f.vB.Kt, data = fishData, prior = pr.vB.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.Kt"))
  
  m.vB.L <- brm(f.vB.L, data = fishData, prior = pr.vB.L, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.L"))
  
  
  m.vB.0
  m.vB <- loo::nlist(m.vB.0, m.vB.K, m.vB.t, m.vB.Kt)
  loo.vB <- m.vB %>% map(loo, cores = 10)
  loo.vB %>% loo_compare
}

#------------

#-----------------
# Gomperts models
#-----------------
if(F){
  f.G.0 <-  bf(sl~ (150 + Linf) * exp(-exp(-(1.28 + K) * (age - t0))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
  f.G.K <-  bf(sl~ (150 + Linf) * exp(-exp(-(1.28 + K) * (age - t0))), K ~ sex, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
  f.G.t <-  bf(sl~ (150 + Linf) * exp(-exp(-(1.28 + K) * (age - t0))), K ~ 1, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
  f.G.Kt <-  bf(sl~ (150 + Linf) * exp(-exp(-(1.28 + K) * (age - t0))), K ~ sex, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
  
  pr.G <- 
    prior(normal(0, 1), nlpar = "K") +  
    prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
    prior(normal(0, 1), nlpar = "t0")
  
  pr.G.K <- pr.G + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") + 
    prior(normal(0, 3), nlpar = "Linf")
  pr.G.t <- pr.G + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") + 
    prior(normal(0, 3), nlpar = "Linf")
  pr.G.Kt <- pr.G + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") +  
    prior(normal(0, 0.2), nlpar = "K", coef = "sex2") + 
    prior(normal(0, 1), nlpar = "Linf")
  
  m.G.0 <- brm(f.G.0, data = fishData, prior = pr.G, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.0"))
  m.G.K <- brm(f.G.K, data = fishData, prior = pr.G.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.K"))
  m.G.t <- brm(f.G.t, data = fishData, prior = pr.G.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.t"))
  m.G.Kt <- brm(f.G.Kt, data = fishData, prior = pr.G.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.Kt"))
  
  m.G <- loo::nlist(m.G.0, m.G.K, m.G.t, m.G.Kt)
  loo.G <- m.G %>% map(loo, cores = 10)
  loo.G %>% loo_compare
}

#------------

#-----------------
# logistic models
#-----------------
if(F){
  f.log.0 <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
  f.log.K <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
  f.log.t <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
  f.log.Kt <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
  
  f.log.L <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ 1, Linf ~ sex + (1|haul), nl = TRUE)
  f.log.LK <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ 1, Linf ~ sex + (1|haul), nl = TRUE)
  f.log.Lt <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ sex, Linf ~ sex + (1|haul), nl = TRUE)
  f.log.LKt <-  bf(sl ~ (152 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ sex, Linf ~ sex + (1|haul), nl = TRUE)
  
  pr.log <- 
    prior(normal(0, 1), nlpar = "K") +  
    prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
    prior(normal(0, 1), nlpar = "t0") +
    prior(normal(0, 10), nlpar = "Linf")
  
  pr.log.K <- pr.log + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") 
  pr.log.t <- pr.log + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
  pr.log.Kt <- pr.log.K + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
  
  pr.log.L <- pr.log + prior(normal(0, 0.2), nlpar = "Linf", coef = "sex2") 
  pr.log.LK <- pr.log.L + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") 
  pr.log.Lt <- pr.log.L + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
  pr.log.LKt <- pr.log.L + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
  
  m.log.0 <- brm(f.log.0, data = fishData, prior = pr.log, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.0"))
  m.log.K <- brm(f.log.K, data = fishData, prior = pr.log.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.K"))
  m.log.t <- brm(f.log.t, data = fishData, prior = pr.log.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.t"))
  m.log.Kt <- brm(f.log.Kt, data = fishData, prior = pr.log.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.Kt"))
  
  m.log.L <- brm(f.log.L, data = fishData, prior = pr.log.L, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.L"))
  m.log.LK <- brm(f.log.LK, data = fishData, prior = pr.log.LK, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.LK"))
  m.log.Lt <- brm(f.log.Lt, data = fishData, prior = pr.log.Lt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.Lt"))
  m.log.LKt <- brm(f.log.LKt, data = fishData, prior = pr.log.LKt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.LKt"))

  m.log.0
  
  m.log <- list()
  for(m in c("0","K","t","Kt","L","LK","Lt","LKt")){
    #m.log[[m]] <- paste0(sv_dir,"/m.log.",m,".rds") %>% readRDS()
  }
  
  loo.log <- m.log %>% map(loo, cores = 10)
  loo.log %>% loo_compare
  
  # combine all models
  loo.m <- c(loo.vB,loo.G,loo.log)
  #saveRDS(loo.m, paste0(sv_dir,"/loo.m.rds"))
}

#-----------------------------------
# Model comparison
#-----------------------------------

loo.m <- readRDS(paste0(sv_dir,"/loo.m.rds"))

make_plot_data <- function(metric_data, levels = names(metric_data)){
  best_model <- metric_data %>% map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% map_dbl(mean),
         metric_diff = metric - max(metric),
         se = metric_data %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_diff = metric_data %>% map(~ .x - metric_data[[best_model]]) %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}

loo_pointwise <- loo.m %>% map("pointwise") %>% map_dfc(~ .[,"elpd_loo"])

se_type <- "se_diff"

loo_df <- make_plot_data(loo_pointwise) %>% 
  mutate(fun = c(rep("vB",4),rep("G",4),rep("log",4)), 
         fe = factor(rep(c("0","K","t","Kt"),3),levels = c("0","K","t","Kt")),
         dim = rep(c(0,1,1,2),3),
         se_ose = se,
         se = .[[se_type]])

loo_ose_models <- loo_df %>% filter(metric + se >= max(metric)) %>% filter(dim == min(dim))
fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")

# faceted plot
loo_plot <- loo_df %>% 
  ggplot(aes(x = fe)) +
  #geom_line(aes(y = metric_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = metric_diff, col = fun), size = 2, show.legend = F) +
  geom_linerange(aes(ymin = metric_diff - se, ymax = metric_diff + se, col = fun), show.legend = F) +
  theme_classic() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90"),
        axis.title.y = element_text(size = 8)) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = metric_diff), shape = 1, size = 6, col = "black", data = loo_ose_models) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  scale_y_continuous(labels=function(x)x*1000, limits = c(-0.022,0.002)) +
  labs(x = NULL, subtitle = "Pinfish PSIS-LOO estimates", y =  expression(Delta*"ELPD "~group("(", 10^- 3,")")))

loo_plot
#ggsave("plots/pinfish_loo_se_mod_2021_10_01.pdf", width = 90, height = 80, units = "mm", device = cairo_pdf()); dev.off()  


#---------------------
# Leave-one-group-out
#---------------------


# fit logo models (i.e., leave-one-haul-out CV)
if(F){
  plan(multisession(workers = 20))
  options(brms.file_refit = "on_change", auto_write = TRUE)
  logo.m <- list()
  for(m in names(loo.m)){
    warning(paste("Fitting model:"), m)
    logo.m[[m]] <- 
      readRDS(paste0(sv_dir,"/",m,".rds")) %>% 
      kfold(chains = 2, future = T, allow_new_levels = T, sample_new_levels = "gaussian", group = "haul")
  }
  #logo.m %>% saveRDS(paste0(sv_dir,"/","logo.m_gaussian",".rds"))
}

if(F){
  plan(multisession(workers = 20))
  options(brms.file_refit = "on_change", auto_write = TRUE)
  logo.m.log <- lapply(m.log, function(m) kfold(m, chains = 2, future = T, allow_new_levels = T, sample_new_levels = "gaussian", group = "haul") )
  logo.m.log %>% saveRDS(paste0(sv_dir,"/","logo.m.log_gaussian",".rds"))
}

m.log$Kt

logo.m.log %>% loo_compare()
logo.m %>% loo_compare()

logo.m <- readRDS(paste0(sv_dir,"/logo.m_gaussian.rds"))

# se_type <- "se_mod" # set above for loo
logo_pointwise <- logo.m %>% map("pointwise") %>% map_dfc(~ .[,"elpd_kfold"])

logo_df <- make_plot_data(logo_pointwise) %>% 
  mutate(fun = c(rep("vB",4),rep("G",4),rep("log",4)), 
         fe = factor(rep(c("0","K","t","Kt"),3),levels = c("0","K","t","Kt")),
         dim = rep(c(0,1,1,2),3),
         se_ose = se,
         se = .[[se_type]])
  
logo_ose_models <- logo_df %>% filter(metric + se >= max(metric)) %>% filter(dim == min(dim))       

fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")


# faceted plot
logo_plot <- logo_df %>% 
  ggplot(aes(x = fe)) +
  #geom_line(aes(y = metric_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = metric_diff, col = fun), size = 2, show.legend = F) +
  geom_linerange(aes(ymin = metric_diff - se, ymax = metric_diff + se, col = fun), show.legend = F) +
  theme_classic() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90"),
        axis.title.y = element_text(size = 8)) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = metric_diff), shape = 1, size = 6, col = "black", data = logo_ose_models) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  scale_y_continuous(labels=function(x)x*1000, limits = c(-0.022,0.002)) +
  labs(x = NULL, subtitle = "Pinfish LOGO-CV estimates", y = expression(Delta*"ELPD "~(10^{"-"~3})))

logo_plot
#ggsave(paste0("plots/pinfish_logo_",se_type,"_2021_10_05.pdf"), width = 90, height = 80, units = "mm", device = cairo_pdf()); dev.off()  


ggpubr::ggarrange(loo_plot,logo_plot, labels = "AUTO", nrow = 1)
#ggsave(paste0("plots/pinfish_ose_",se_type,"_2021_10_05.pdf"), width = 180, height = 80, units = "mm", device = cairo_pdf()); dev.off()  



#--------------
## DATA PLOTS
#--------------


# Plot 1 - conditional model

m.vB.0 <- readRDS(paste0(sv_dir,"/m.vB.0.rds"))
hauls <- m.vB.0$data$haul %>% unique
ce_data <- hauls %>% map_dfr(~ tibble(age = seq(0.5,6.3,0.04), haul = .x))

g.vB <- function(age, Linf, K, t0) {(150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0))))}

#pp_m.vB.0 <- m.vB.0 %>% posterior_predict(ndraws = 500)

ps.0 <- m.vB.0 %>% as_tibble %>% 
  rename_with(~ str_replace(.x, fixed("r_haul__Linf["),""), starts_with("r_haul")) %>% 
  rename_with(~ str_replace(.x, fixed(",Intercept]"),""), ends_with("Intercept]")) %>% 
  rename(K = b_K_Intercept, t0 = b_t0_Intercept, Linf0 = b_Linf_Intercept, tau = sd_haul__Linf_Intercept) %>% 
  select(-lp__) %>% sample_n(500) %>% mutate(s = 1:500)

ce_plot_data <- ps.0 %>% pivot_longer(-any_of(c("sigma","tau","s","K","Linf0","t0")),names_to = "haul", values_to = "Linf_haul") %>% 
  mutate(Linf = Linf0 + Linf_haul) %>% 
  full_join(ce_data, by = "haul") %>% 
  group_by(s) %>% 
  mutate(sig_err = rnorm(1,0,mean(sigma))) %>%  ## generate residual error
  mutate(length = g.vB(age,Linf,K,t0) + sig_err) 
  
ribbon_data <- ce_plot_data %>% 
  group_by(age) %>% 
  summarise(l_low = quantile(length, 0.025), l_hi = quantile(length, 0.975)) 
#  group_by(age) %>% 
#  summarise(l_low = min(l_low), l_hi = max(l_hi))

ce_plot_data_col <- ce_plot_data %>% 
  group_by(haul) %>% 
  mutate(haul_col = mean(Linf) + 150) 

ce_plot_data_col
obs_data <- m.vB.0$data %>% left_join(ce_plot_data_col %>% select(haul, haul_col) %>% distinct, by = "haul")

ce_plot_cond <- ce_plot_data_col %>% 
  group_by(haul_col,age) %>% 
  summarise(length = mean(length)) %>% 
  ggplot(aes(x = age)) +
  #geom_ribbon(aes(ymin = l_low, ymax = l_hi), fill = "black", alpha = 0.1, data = ribbon_data) +
  geom_line(aes(y = length, group = haul_col, col = haul_col), size = 0.25, alpha = 0.4) +
  geom_point(aes(y = sl, group = haul_col, col = haul_col), size = 0.5, alpha = 1, data = obs_data) +
  #geom_point(aes(y = sl), size = 0.3, alpha = 1, col = "white", data = m.vB.0$data) +
  theme_classic() +
  theme(legend.position = "bottom", legend.key.height =  unit(5,"mm")) +
  labs(col = expression(L[0]~"| haul"~phantom(0))) +
  guides(color=guide_colourbar(title.vjust=0.5))

ce_plot_cond

# Plot 2 - marginal model
m.log.t <- readRDS(paste0(sv_dir,"/m.log.t.rds"))
hauls <- m.log.t$data$haul %>% unique
ce_data_2 <- c(1,2) %>% map_dfr(~ tibble(age = seq(0.5,6.3,0.04), sex = .x))

#m.log.t$formula
g.log.t <- function(age, Linf, K, t0, t0sex) {(151 + Linf)/(1 + exp(-(1.4 + K) * (age - (t0 + t0sex))))}

ps.t <- m.log.t %>% as_tibble %>% 
  select(K = b_K_Intercept, t0 = b_t0_Intercept, Linf0 = b_Linf_Intercept, 
         tau = sd_haul__Linf_Intercept, t0_sex2 = b_t0_sex2, sigma) %>% 
  sample_n(500) %>% mutate(s = 1:500)

ce_plot_data2 <- c(1,2) %>% map_dfr(~ ps.t %>% mutate(sex = .x)) %>% 
  full_join(ce_data_2, by = "sex") %>% 
  group_by(s) %>% 
  mutate(sig_err = rnorm(1,0,mean(sigma))) %>%  ## draw from residual error model
  mutate(tau_err = rnorm(1,0,mean(tau))) %>% ## draw from hierarchical model for Linf
  mutate(t0_sex2 = t0_sex2*(sex - 1)) %>% 
  mutate(length_marg = g.log.t(age,Linf0 + tau_err,K,t0,t0_sex2) + sig_err )
  
# select colours
dark2 <- RColorBrewer::brewer.pal(8,"Dark2")
pie(rep(1, length(dark2)), col = dark2 , main="") 
  
ce_plot_marg <- ce_plot_data2 %>% 
  mutate(length = length_marg) %>% 
  group_by(age,sex) %>% 
  summarise(l_low = quantile(length, 0.025),
            l_hi = quantile(length, 0.975),
            length = mean(length)) %>% 
  mutate(sex = factor(sex)) %>% 
  ggplot(aes(x = age)) +
  geom_ribbon(aes(ymin = l_low, ymax = l_hi), fill = blues9[3], alpha = 0.4) +
  geom_line(aes(y = length, col = sex), size = 0.5, lty = "solid") +
  geom_point(aes(y = sl, col = sex), size = 0.4, alpha = 1, data = m.log.t$data) +
  scale_color_manual(values = c(dark2[1],dark2[8])) +
  #geom_point(aes(y = sl), size = 0.3, alpha = 1, col = "white", data = m.vB.0$data) +
  theme_classic() +
  theme(legend.position = "bottom", legend.key.height =  unit(10,"mm"))

yl2 <- ylim(c(45,206))

ggpubr::ggarrange(ce_plot_cond + labs(subtitle = "Pinfish data: conditional focus") + yl2, 
                  ce_plot_marg + yl2 + labs(y = "", subtitle = "Pinfish data: marginal focus"), nrow = 1, labels = "AUTO")

ggsave(paste0("plots/pinfish_data_2_2021_10_06.pdf"), width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()  


ggpubr::ggarrange(loo_plot + labs(subtitle = "PSIS-LOO-CV estimates"),
                  logo_plot + labs(subtitle = "LOGO-CV estimates", y = " "),
                  ce_plot_cond + labs(subtitle = "Conditional prediction: vB|0") + yl2, 
                  ce_plot_marg + yl2 + labs(y = "", subtitle = "Marginal prediction: log|t"), 
                  nrow = 2, ncol = 2, labels = "AUTO", heights = c(2,3))

ggsave(paste0("plots/pinfish_panel_2021_10_06.png"), width = 180, height = 180, units = "mm", device = cairo_pdf()); dev.off()  

