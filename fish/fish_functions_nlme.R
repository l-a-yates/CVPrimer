#------------------------------------------------------
# Functions for pinfish non-linear mixed effects models
# Luke Yates
#------------------------------------------------------

library(tidyverse)
library(fishmethods)
library(nlme)
library(pbmcapply)
library(knitr)

MAX_CORES <- 20

#---------------------------------
# Part 1: Haul as a random effect
#---------------------------------

# fit full set of random effects models
fit.re <- function(data){
  
  fits <- list()

  # von Bertalanffy 
  fits$vB.0 <-  nlme(sl ~ Linf * (1 - exp(-(K * (age - t0)))),
       fixed = list(Linf + K + t0 ~ 1),
       random = Linf ~ 1|haul,
       data = data,
       start = c(Linf = 210, K = 0.3, t0 = -1))
  
  fits$vB.K <-  nlme(sl ~ Linf * (1 - exp(-(K * (age - t0)))),
                fixed = list(K ~ sex, Linf + t0 ~ 1),
                random = Linf ~ 1|haul,
                data = data,
                start = c(K1 = 0.3, K2 = 0.3, Linf = 210, t0 = -1))
  
  fits$vB.t <-  nlme(sl ~ Linf * (1 - exp(-(K * (age - t0)))),
                fixed = list(t0 ~ sex, Linf + K ~ 1),
                random = Linf ~ 1|haul,
                data = data,
                start = c(t01 = -1, t02 = -1,Linf = 210, K = 0.3))
  

  fits$vB.Kt <- nlme(sl ~ Linf * (1 - exp(-(K * (age - t0)))),
                fixed = list(K + t0 ~ sex, Linf ~ 1),
                random = Linf ~ 1|haul,
                data = data,
                start = c(K1 = 0.3, K2 = 0.3, t01 = -1, t02 = -1, Linf = 210))
  
  # Gomperts
  fits$G.0 <-  nlme(sl~ Linf * exp(-exp(-K * (age - t0))),
                     fixed = list(Linf + K + t0 ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(Linf = 210, K = 0.3, t0 = -1))
  
  fits$G.K <-  nlme(sl~ Linf * exp(-exp(-K * (age - t0))),
                     fixed = list(K ~ sex, Linf + t0 ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(K1 = 0.3, K2 = 0.3, Linf = 210, t0 = -1))
  
  fits$G.t <-  nlme(sl~ Linf * exp(-exp(-K * (age - t0))),
                     fixed = list(t0 ~ sex, Linf + K ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(t01 = -1, t02 = -1,Linf = 210, K = 0.3))
  
  
  fits$G.Kt <-  nlme(sl~ Linf * exp(-exp(-K * (age - t0))),
                      fixed = list(K + t0 ~ sex, Linf ~ 1),
                      random = Linf ~ 1|haul,
                      data = data,
                      start = c(K1 = 0.3, K2 = 0.3, t01 = -1, t02 = -1, Linf = 210))
  
  # logistic
  fits$log.0 <-  nlme(sl ~ Linf/(1 + exp(-K * (age - t0))),
                     fixed = list(Linf + K + t0 ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(Linf = 210, K = 0.3, t0 = -1))
  
  fits$log.K <-  nlme(sl ~ Linf/(1 + exp(-K * (age - t0))),
                     fixed = list(K ~ sex, Linf + t0 ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(K1 = 0.3, K2 = 0.3, Linf = 210, t0 = -1))
  
  fits$log.t <-  nlme(sl ~ Linf/(1 + exp(-K * (age - t0))),
                     fixed = list(t0 ~ sex, Linf + K ~ 1),
                     random = Linf ~ 1|haul,
                     data = data,
                     start = c(t01 = -1, t02 = -1,Linf = 210, K = 0.3))
  
  
  fits$log.Kt <-  nlme(sl ~ Linf/(1 + exp(-K * (age - t0))),
                      fixed = list(K + t0 ~ sex, Linf ~ 1),
                      random = Linf ~ 1|haul,
                      data = data,
                      start = c(K1 = 0.3, K2 = 0.3, t01 = -1, t02 = -1, Linf = 210))
  
  fits
}

# compute point-wise conditional log-likelihood for nlme and nls fits
comp_ll <- function(fit, newdata = NULL){
  sigma <- sqrt(mean(residuals(fit)^2))
  if(is.null(newdata)){
    mu = residuals(fit)
    y = 0
  } else {
    mu = predict(fit, newdata = newdata)
    y = newdata %>% pull('sl')
  }
  dnorm(y,mu,sigma,log = T) %>% as.numeric
}

# compute loo and summary statistics for a supplied data set and fitting function
loocv <- function(data, fit_fun, bias.adjust = T){
  ld_i_ws = map_dfc(fit_fun(data), comp_ll, newdata = data)
  cv_i = pbmclapply(1:nrow(data), function(i){
    fit.train = fit_fun(data[-i,])
    adj = sapply(fit.train, function(fit) comp_ll(fit, newdata = data) %>% mean)
    ld = sapply(fit.train, comp_ll, newdata = data[i,])
    list(adj = adj, ld = ld)
  }, mc.cores = MAX_CORES)
  adj_i = cv_i %>% map_df("adj")
  bias_i = (ld_i_ws - adj_i)
  ld_i = cv_i %>% map_df("ld")
  if(bias.adjust) ld_i = ld_i + bias_i
  eld = ld_i %>% map_dbl(mean)
  se = ld_i %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  best_model = ld_i %>% map_dbl(mean) %>% which.max() %>% names()
  eld_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(mean)
  se_mod = sqrt(1 -cor(ld_i)[best_model,])*se[best_model] 
  se_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  p_loo = map_dbl(ld_i_ws - ld_i, sum); p_loo
  tibble(model = names(eld), eld_diff = eld_diff, se_mod = se_mod,
         se_diff = se_diff, eld = eld, se = se, p_loo = p_loo)
}

# compute K-fold CV and summary statistics for a supplied data set and fitting function
kfcv <- function (data, fit_fun, K = 10, bias.adjust = T, maxIter = 5){
  # fold data and check that all training fits converge
    if(K>10) stop("K exceeds minimum group size. Try K <= 10.")
    fitted <- F; it <- 0
    while(!fitted){
      message("Folding data and checking convergence of model fits to training data.")
      f_data <- group_by(data,haul) %>% mutate(fold = sample(rep(sample(K), length = n()))) %>% 
        ungroup() %>% arrange(fold)
      it <- it + 1
      try(fits <- pbmclapply(1:K, function(k) fit.re(filter(f_data,fold!=k)), mc.cores = min(K,MAX_CORES)))
      fitted <- exists("fits")
      if (it > maxIter) stop("MaxIter exceeded")
    }
  ld_i_ws = map_dfc(fit_fun(f_data), comp_ll, newdata = f_data) # within-sample fit
  cv_i = lapply(1:K, function(k){ # point-wise cv fits
    ld_adj = sapply(fits[[k]], comp_ll, newdata = f_data) %>% as.data.frame %>% mutate(obs = 1:nrow(data))
    ld = sapply(fits[[k]], comp_ll, newdata = filter(f_data,fold==k)) %>% as.data.frame
    list(ld_adj = ld_adj, ld = ld)
  })
  ld_i = cv_i %>% map_dfr("ld")
  ld_adj_i = cv_i %>% map_dfr("ld_adj") %>% group_by(obs) %>% summarise(across(everything(),mean)) %>% select(-obs)
  bias_i = (ld_i_ws - ld_adj_i)
  if(bias.adjust) ld_i = ld_i + bias_i
  eld = ld_i %>% map_dbl(mean)
  se = ld_i %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  best_model = eld %>% which.max() %>% names()
  eld_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(mean)
  se_mod = sqrt(1 -cor(ld_i)[best_model,])*se[best_model] 
  se_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  p_loo = map_dbl(ld_i_ws - ld_i, sum)
  tibble(model = names(eld), eld_diff = eld_diff, se_mod = se_mod,
         se_diff = se_diff, eld = eld, se = se, p_loo = p_loo)  
}


#------------------------------------------
# Part 2: Haul as a random vs fixed effects
#------------------------------------------

# simple vB fixed effects model, ignore sex
fit.fe.0 <- function(data){
  n.haul <- length(unique(data$haul))
  vb.fe <- sl~Linf[haul]*(1-exp(-K*(age-t0)))
  sv.fe <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,1,1))
  nls(vb.fe,data=data,start=sv.fe)
}


# simple vB random effects model, ignore sex
fit.re.0 <- function(data){
  nlme(sl ~ Linf * (1 - exp(-(K * (age - t0)))),
       fixed = list(Linf + K + t0 ~ 1),
       random = Linf ~ 1|haul,
       data = data,
       start = c(Linf = 210, K = 0.3, t0 = -1))
}


# fit the fe and re models
fit.re.fe.0 <- function(data) list(fe = fit.fe.0(data), re = fit.re.0(data))
