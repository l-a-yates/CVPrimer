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
#library(FlexParamCurve)
library(RColorBrewer)
library(brms)

rm(list=ls())

#source("fish_functions_nlme.R")
options(brms.file_refit = "on_change", auto_write = TRUE)

run_date <- "2021_10_07"
#sv_dir <- paste0("fits_",run_date)
#if(!dir.exists(sv_dir)) dir.create(sv_dir)
#sv_dir <- paste0("fits_",run_date)

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
min.size <- 5 # set minimum haaul size
fishData <- fishData %>% group_by(haul) %>% filter(n()>= min.size) %>% ungroup
fishData <- fishData %>% filter(haul != "TBD930066") # remove outlier
fishData$haul <-  factor(fishData$haul)


#----------------------
# Define and fit models
#----------------------

# growth functions
nlf.vB <- sl ~ (152 + Linf) * (1 - exp(-1*((1.3 + K) * (age - t0))))
nlf.G <- sl~ (152 + Linf) * exp(-exp(-(1.3 + K) * (age - t0)))
nlf.log <- sl ~ (152 + Linf)/(1 + exp(-(1.3 + K) * (age - t0)))

# specify brms formula
make_form <- function(fun, K, t, L) {
  if(K) {f.K <- K ~ sex} else {f.K <- K ~ 1}
  if(t) {f.t <- t0 ~ sex} else {f.t <- t0 ~ 1}
  if(L) {f.L <- Linf ~ sex + (1|haul)} else {f.L <- Linf ~ 1 + (1|haul)}
  bf(get(paste0("nlf.",fun)), f.K, f.t, f.L, nl = T)
}

# specify priors
make_prior <- function(K, t, L, sig_Linf = 10, sig_sex2 = 0.3){
  p <-  
    prior(normal(0, 1), nlpar = "K") +  
    prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
    prior(normal(0, 1), nlpar = "t0") +
    set_prior(paste0("normal(0, ",sig_Linf,")"), nlpar = "Linf")
  if(K) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "K", coef = "sex2") 
  if(t) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "t0", coef = "sex2") 
  if(L) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "Linf", coef = "sex2") 
  return(p)
}

# characterise model set
models_grid <- expand.grid(fun = c("vB","G","log"), K = c(0,1), L = c(0,1), t  = c(0,1), stringsAsFactors = F) %>% 
  mutate(name = paste0(fun,".",ifelse(K,"K",""),ifelse(L,"L",""),ifelse(t,"t",""), ifelse(K + L + t, "","0")))

# fit and save models
if(F){
  seed <- 60869 #sample(1e5,1)
  m <- list()
  plan(multisession(workers = 24))
  for(i in 1:nrow(models_grid)){
    print(paste("Fitting model", i, models_grid[i,"name"])) 
    m[[models_grid[i,"name"]]] <- futureCall(brm, args = list(formula = with(models_grid[i,], make_form(fun,K,t,L)),
                                                              prior = with(models_grid[i,], make_prior(K,t,L)),
                                                              data = fishData,
                                                              future = F,
                                                              chains = 4,
                                                              seed = seed,
                                                              iter = 4000, 
                                                              file = paste0(sv_dir,"/",models_grid[i,"name"]),
                                                              file_refit = "always"),
                                             seed = T,
                                             earlySignal = T
    )
  }
}

# load all models and check convergence
m.fits <- models_grid %>% pull(name, name = name) %>% map(~ paste0(sv_dir,"/",.x,".rds") %>% readRDS)
m.fits %>% map_dbl(~ .x %>% rhat() %>% map_dbl(mean) %>% mean)

# 2 models have convergence issues: refit them with narrower priors
if(F){
  m.fits$vB.KLt %>% update(cores = 4, inits = "0", seed = seed, prior = make_prior(1,1,1,sig_Linf = 5, sig_sex2 = 0.2), file = paste0(sv_dir,"/vB.KLt_update"))
  m.fits$vB.Kt %>% update(cores = 4, inits = "0", seed = seed, prior = make_prior(1,1,0,sig_Linf = 5, sig_sex2 = 0.2), file = paste0(sv_dir,"/vB.Kt_update"))
}

# refits successful: replace them in the list of fits
m.fits$vB.KLt <- readRDS(paste0(sv_dir,"/vB.KLt_update.rds"))
m.fits$vB.Kt <- readRDS(paste0(sv_dir,"/vB.Kt_update.rds"))

# compute PSIS-LOO
if(F){
  m.loo <- m.fits %>% furrr::future_map(loo)
  #saveRDS(m.loo,paste0(sv_dir,"/m.loo.rds"))
}
m.loo <- readRDS(paste0(sv_dir,"/m.loo.rds"))

# compute LOGO CV estimates
if(F){
  plan(multisession(workers = 45))
  m.logo.cv <- lapply(m.fits, kfold, chains = 2, future = T, group = "haul")
  #m.logo.cv %>% saveRDS(paste0(sv_dir,"/m.logo.rds"))
}

m.logo.cv<- readRDS(paste0(sv_dir,"/m.logo.rds"))


#------
# Plots
#------

# extract pointwise loo estimates
model_names <- models_grid %>% mutate(dim = K + L + t) %>% arrange(fun,dim) %>% pull(name)
m.loo.pointwise <- m.loo %>% map("pointwise") %>% map_dfc(~ .x[,"elpd_loo"]) %>% relocate(all_of(model_names))

m.best <- m.loo.pointwise %>% map_dbl(mean) %>% which.max() %>% names
m.loo.pointwise %>% mutate(across(everything(), ~ .x - m.loo.pointwise[[m.best]])) %>% 
  {tibble(model = names(.) %>% {factor(., levels = .)},
          elpd_diff = map_dbl(.,mean),
          se_diff = map_dbl(.,sd)/sqrt(nrow(.)))} %>% 
ggplot(aes(model)) +
  geom_point(aes(y =  elpd_diff))+
  geom_linerange(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff)) +
  theme_classic()



# LOGO
m.logo.pointwise <- m.logo.cv %>% map("pointwise") %>% map_dfc(~ .x[,"elpd_kfold"]) %>% relocate(all_of(model_names))
m.logo.best <- m.logo.pointwise %>% map_dbl(mean) %>% which.max() %>% names

m.logo.pointwise %>% mutate(across(everything(), ~ .x - m.logo.pointwise[[m.logo.best]])) %>% 
  {tibble(model = names(.) %>% {factor(., levels = .)},
          elpd_diff = map_dbl(.,mean),
          se_diff = map_dbl(.,sd)/sqrt(nrow(.)))} %>% 
  filter(elpd_diff > -0.03) %>% 
  ggplot(aes(model)) +
  geom_point(aes(y =  elpd_diff))+
  geom_linerange(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff)) +
  theme_classic()

