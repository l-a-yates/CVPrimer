#------------------------------------------------------
# Functions for pinfish non-linear mixed effects models
# Luke Yates
# last update 18/02/2020
#------------------------------------------------------

library(tidyverse)
library(fishmethods)
library(nlme)
library(pbmcapply)
library(knitr)


# von Bertalanffy growth function
#vb <- function(age, Linf, K, t0) Linf * (1 - exp(-(K * (age - t0))))

MAX_CORES <- 20

#-----------------------------
# Random Effects Only - part 1
#-----------------------------


# fit set of random effects models
fit.re <- function(data){
  
  fits <- list()

  # vB 
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
  
  #logisitc
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

# compute loo and summary statistics for a supplied data set, fitting function, and valid score type
loocv <- function(data, fit_fun, bias.adjust = T, score = "eld"){
  if(score == "eld"){
    loocv.eld(data = data, fit_fun = fit_fun, bias.adjust = bias.adjust)
  } else if(score == "mse"){
    loocv.mse(data = data, fit_fun = fit_fun, bias.adjust = bias.adjust)
  } else warning("Select a valid score.")
}



# compute loo and summary statistics for a supplied data set and fitting function
loocv.mse <- function(data, fit_fun, bias.adjust = T){
  ws_rs = map_dfc(fit_fun(data), ~ .x %>% residuals %>% as.numeric %>% {.^2})
  cv_i = pbmclapply(1:nrow(data), function(i){
    fit.train <- fit_fun(data[-i,])
    adj = sapply(fit.train, function(fit) sum((as.numeric(predict(fit, newdata = data)) - data[,'sl'])^2)/nrow(data))
    rs = sapply(fit.train, function(fit) {(as.numeric(predict(fit, newdata = data[i,]) - data[i,'sl']))^2})
    list(rs = rs, adj = adj)
  }, mc.cores = MAX_CORES)
  adj_i = cv_i %>% map_df("adj") 
  bias_i = (ws_rs - adj_i)
  rs_i <- cv_i %>% map_df("rs")
  if(bias.adjust) rs_i <- rs_i + bias_i
  mse <- rs_i %>% map_dbl(mean)
  se_ose <- rs_i %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  best_model <- rs_i %>% map_dbl(mean) %>% which.min() %>% names()
  delta_mse <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(mean)
  se_mod <- sqrt(1 -cor(rs_i)[best_model,])*se_ose[best_model] 
  se_diff <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  p_loo <- log(mse/map_dbl(ws_rs,mean))*(nrow(data)/2)
  tibble(model = names(mse), delta_mse = delta_mse, se_mod = se_mod,
         se_diff = se_diff, mse = mse, se_ose = se_ose, p_loo = p_loo)
}


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
loocv.eld <- function(data, fit_fun, bias.adjust = T){
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


kfcv <- function (data, fit_fun, K = 10, bias.adjust = T, maxIter = 5){
  # fold data and check that all training fits converge
    if(K>10) stop("K exceeds group size. Try K <= 10.")
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
  ld_i_ws = map_dfc(fit_fun(f_data), comp_ll, newdata = f_data)
  cv_i = lapply(1:K, function(k){
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
  best_model = ld_i %>% map_dbl(mean) %>% which.max() %>% names()
  eld_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(mean)
  se_mod = sqrt(1 -cor(ld_i)[best_model,])*se[best_model] 
  se_diff = ld_i %>% map(~ .x - ld_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  p_loo = map_dbl(ld_i_ws - ld_i, sum)
  tibble(model = names(eld), eld_diff = eld_diff, se_mod = se_mod,
         se_diff = se_diff, eld = eld, se = se, p_loo = p_loo)  
}

# Implements K-fold cross validation
# Bias adjustment follows Davison and Hinkley (Due originally to Burman 1989)
# W.r.t variance, it is better to increase K and decrease REP than vice versa (Burman 1989)
# Function output includes CV-scores (out-of-sample deviance), frequency selection probabilities and 
# complexity estimates using: complexity = ll_ws - ll_oos (i.e. difference between within- and out-of- 
# sample log-likelihood)
# K-Fold performs better than Leave-g-out (g = n/k) as an estimator of out-of-sample deviance (Expected KL discrepancy)
# na.rm = T is not ideal, but help with K = 2 predictions which occasionally fail (~1% of iterations)
re.kfcv <- function (data, K = 10, REP = 1, maxIter = 5, bias.adj = T, na.rm = F){
  n <- nrow(data)
  ws.mc <- sapply(fit.re(data), function(m) mean(as.numeric(m$residuals[,2])^2))
  CV.iter <- pbmclapply(1:REP, function(r){
    
    fitted <- F; it <- 0
    while(!fitted){
      fold <- sample(rep(1:K, length = n))
      it <- it + 1
      try(fits <- lapply(1:K, function(k) fit.re(data[fold!=k,])))
      fitted <- exists("fits")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    
    mc <- adj <- vector("list",K)
    
    for(k in 1:K){
      test <- (fold == k)
      mc[[k]] <- sapply(fits[[k]], function(fit) {(as.numeric(predict(fit, newdata = data[test,]) - pull(data[test,'sl'])))^2})
      adj[[k]] <- sapply(fits[[k]], function(fit) sum((as.numeric(predict(fit, newdata = data)) - pull(data[,'sl']))^2)*sum(test)/nrow(data)^2)
    }
    mck.ba <- lapply(1:K, function(i) apply(mc[[i]],2,mean) + ws.mc - (adj[[i]]*n)/sum(fold == i))
    fsel <- apply(sapply(mc, function(m) apply(m,2,mean) == min(apply(m,2,mean))),1,sum)/K
    fsel.ba <- apply(sapply(mck.ba, function(x) x == min(x)),1,sum)/K
    mci <- bind_rows(lapply(mc, as.data.frame))
    #return(list(mc=mc, adj=adj))
    mc <- apply(mci, 2 , mean)
    err <- apply(mci, 2, sd)/sqrt(n)
    if(bias.adj) mc <- mc + ws.mc -apply(sapply(adj, function(x) x),1,sum)
    list(mc = mc, err = err, fsel = fsel, fsel.ba = fsel.ba)
  }, mc.cores = min(REP,MAX_CORES))
  #return(CV.iter)
  mc <- apply(sapply(CV.iter, function(x) x$mc),1,mean, na.rm = na.rm)
  err <- sqrt(apply(sapply(CV.iter, function(x) (x$err)^2),1,mean, na.rm = na.rm))
  fsel <-  apply(sapply(CV.iter, function(x) x$fsel),1,sum, na.rm = na.rm)
  fsel.ba <-  apply(sapply(CV.iter, function(x) x$fsel.ba),1,sum, na.rm = na.rm)
  complexity <- log(mc/ws.mc)*(nrow(data)/2)
  tibble(Model = names(mc), mc = mc, err = err, Score = n*log(mc), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel), 
         fsel.ba = fsel.ba/sum(fsel.ba), complexity = complexity)
} # end re.kfcv










#--------------------------------
# Random Effects vs Fixed Effects
#--------------------------------



# simple vB fixed effects model, ignore sex
fit.fe.0 <- function(data){
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


re.loocv2 <- function(data, fit_fun, bias.adjust = T){
  ws_rs = map_dfc(fit_fun(data), ~ .x %>% residuals %>% as.numeric %>% {.^2})
  cv_i = pbmclapply(1:nrow(data), function(i){
    fit.train <- fit_fun(data[-i,])
    adj = sapply(fit.train, function(fit) sum((as.numeric(predict(fit, newdata = data)) - data[,'sl'])^2)/nrow(data))
    rs = sapply(fit.train, function(fit) {(as.numeric(predict(fit, newdata = data[i,]) - data[i,'sl']))^2})
    list(rs = rs, adj = adj)
  }, mc.cores = 30)
  adj_i = cv_i %>% map_df("adj") 
  bias_i = (ws_rs - adj_i)
  rs_i <- cv_i %>% map_df("rs")
  if(bias.adjust) rs_i <- rs_i + bias_i
  mse <- rs_i %>% map_dbl(mean)
  se_ose <- rs_i %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  best_model <- rs_i %>% map_dbl(mean) %>% which.min() %>% names()
  delta_mse <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(mean)
  se_mod <- sqrt(1 -cor(rs_i)[best_model,])*se_ose[best_model] 
  se_diff <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
  p_loo <- log(mse/map_dbl(ws_rs,mean))*(nrow(data)/2)
  tibble(model = names(mse), delta_mse = delta_mse, se_mod = se_mod,
         se_diff = se_diff, mse = mse, se_ose = se_ose, p_loo = p_loo)
}

# leave-one-out CV
re.fe.loocv <- function(data){
  CV.iter <- pbmclapply(1:nrow(data), function(i){
    fit.train <- fit.re.fe.0(data[-i,])
    fe <- (as.numeric(predict(fit.train[[1]], newdata = data[i,])) - data[i,'sl'])^2
    re <- (as.numeric(predict(fit.train[[2]], newdata = data[i,])) - data[i,'sl'])^2
    c(fe= fe,re = re)
  }, mc.cores = 30)
  models <- names(CV.iter[[1]])
  CV.iter <- sapply(CV.iter, function(x) as.numeric(x))#; return(CV.iter)
  MSE <- apply(CV.iter,1,mean)
  MSE.se <- apply(CV.iter,1,sd)/sqrt(nrow(data))
  fsel <- apply(apply(CV.iter,2, function(b) b == min(b)),1,sum)
  tibble(Model = models, MSE = MSE, MSE.se = MSE.se, Score = nrow(data)*log(MSE), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel))
}


# leave-g-out CV
re.fe.lgo <- function(g, iter, data, maxIter = 1000,   n.min = 3){
  # p gives sampling probabilities proportionate to haul size (guarantees at least n.min random training samples in each haul)
  p <- data %>% group_by(haul) %>% mutate(p = as.numeric(!1:n() %in% sample(n(),n.min))) %>% pull(p) 
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; it <- 0
    while(!fitted){
      it <- it + 1
      test <- sample(1:nrow(data), g, prob = p)
      try(fit.train <- fit.re.fe.0(data[-test,]))
      fitted <- exists("fit.train")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/g
    re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,])) - data[test,'sl'])^2)/g
    c(fe= fe,re = re)
  }, mc.cores = 50)
  CV.iter <- sapply(CV.iter, function(x) x)#; return(CV.iter)
  MSE <- apply(CV.iter,1,mean)
  MSE.se <- apply(CV.iter,1,sd)*sqrt(g)/sqrt(nrow(data))
  fsel <- apply(apply(CV.iter,2, function(b) b == min(b)),1,sum)
  tibble(Model = rownames(CV.iter), MSE = MSE, MSE.se = MSE.se, Score = nrow(data)*log(MSE), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel))
}


# leave-one-out CV
re.fe.loocv_old <- function(data){
  CV.iter <- pbmclapply(1:nrow(data), function(i){
    fit.train <- fit.re.fe.0(data[-i,])
    fe <- (as.numeric(predict(fit.train[[1]], newdata = data[i,])) - data[i,'sl'])^2
    re <- (as.numeric(predict(fit.train[[2]], newdata = data[i,])) - data[i,'sl'])^2
    c(fe= fe,re = re)
  }, mc.cores = 30)
  models <- names(CV.iter[[1]])
  CV.iter <- sapply(CV.iter, function(x) as.numeric(x))#; return(CV.iter)
  MSE <- apply(CV.iter,1,mean)
  MSE.se <- apply(CV.iter,1,sd)/sqrt(nrow(data))
  fsel <- apply(apply(CV.iter,2, function(b) b == min(b)),1,sum)
  tibble(Model = models, MSE = MSE, MSE.se = MSE.se, Score = nrow(data)*log(MSE), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel))
}





# loocv per group (haul)
re.fe.loopg <- function(iter, data, maxIter = 1000){
  n.group <- length(unique(data$haul))
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; iter <- 0
    while(!fitted){
      iter <- iter + 1
      test <- data %>% group_by(haul) %>% mutate(p = (1:n() %in% sample(n(),1))) %>% pull(p)
      try(fit.train <- fit.re.fe(data[-test,]))
      fitted <- exists("fit.train")
      if (iter > maxIter) stop("MaxIter exceeded")
    }
    fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/n.group
    re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,], level = 1)) - data[test,'sl'])^2)/n.group
    c(fe= fe,re = re)
  }, mc.cores = 50)
  CV.iter <- sapply(CV.iter, function(x) x)#; return(CV.iter)
  MSE <- apply(CV.iter,1,mean)
  MSE.se <- apply(CV.iter,1,sd)*sqrt(n.group)/sqrt(nrow(data))
  fsel <- apply(apply(CV.iter,2, function(b) b == min(b)),1,sum)
  tibble(Model = rownames(CV.iter), MSE = MSE, MSE.se = MSE.se, Score = nrow(data)*log(MSE), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel))
}


# bias-adjusted lgo CV
adj.lgo <- function(g, iter, data, maxIter = 1000){
  n.min <- 3 # min samples left in each haul
  p <- data %>% group_by(haul) %>% mutate(p = as.numeric(!1:n() %in% sample(n(),n.min))) %>% pull(p)   # sampling probabilities 
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; it <- 0
    while(!fitted){
      it <- it + 1
      test <- sample(1:nrow(data), g, prob = p)
      try(fit.train <- fit.re.fe(data[-test,]))
      fitted <- exists("fit.train")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    mse.fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/g
    adj2.fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data)) - data[,'sl'])^2)/nrow(data)
    mse.re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,])) - data[test,'sl'])^2)/g
    adj2.re <- sum((as.numeric(predict(fit.train[[2]], newdata = data)) - data[,'sl'])^2)/nrow(data)
    #re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,], level = 1)) - data[test,'sl'])^2)/g
    list(fe = c(mse.fe = mse.fe, adj2.fe = adj2.fe), re = c(mse.re= mse.re, adj2.re = adj2.re))
  }, mc.cores = 50)
  return(CV.iter)
  CV.iter <- sapply(CV.iter, function(x) x); return(CV.iter)
}






##-----------
## OTHER FUNS
##-----------

if(F){

  
  # leave-one-out CV
  # implements optional bias adjustment by default
  re.loocv_old1 <- function(data, bias.adjust = T){
    ws_rs = sapply(fit.re(data), function(m) as.numeric(m$residuals[,2])^2)
    cv_i = pbmclapply(1:nrow(data), function(i){
      fit.train <- fit.re(data[-i,])
      adj = sapply(fit.train, function(fit) sum((as.numeric(predict(fit, newdata = data)) - data[,'sl'])^2)/nrow(data))
      rs = sapply(fit.train, function(fit) {(as.numeric(predict(fit, newdata = data[i,]) - data[i,'sl']))^2})
      list(rs = rs, adj = adj)
    }, mc.cores = 30)
    adj_i = cv_i %>% map_df("adj") 
    bias_i = (ws_rs - adj_i)
    rs_i <- cv_i %>% map_df("rs")
    if(bias.adjust) rs_i <- rs_i + bias_i
    mse <- rs_i %>% map_dbl(mean)
    se_ose <- rs_i %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
    best_model <- rs_i %>% map_dbl(mean) %>% which.min() %>% names()
    delta_mse <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(mean)
    se_mod <- sqrt(1 -cor(rs_i)[best_model,])*se_ose[best_model] 
    se_diff <- rs_i %>% map(~ .x - rs_i[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x)))
    p_loo <- log(mse/mean(ws_rs))*(nrow(data)/2)
    tibble(model = names(mse), delta_mse = delta_mse, se_mod = se_mod,
           se_diff = se_diff, mse = mse, se_ose = se_ose, p_loo = p_loo)
  }
  
  
# leave-g-out CV, RE
# p gives sampling probabilities proportionate to haul size (guarantees at least n.min random training samples in each haul)
re.lgo <- function(g, iter, data, maxIter = 1000, n.min = 3, bias.adj = T){
  n <- nrow(data)
  ws.mc <- sapply(fit.re(data), function(m) mean(as.numeric(m$residuals[,2])^2))
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; it <- 0
    while(!fitted){
      p <- data %>% group_by(haul) %>% mutate(p = as.numeric(!1:n() %in% sample(n(),n.min))) %>% pull(p) 
      it <- it + 1
      test <- sample(1:nrow(data), g, prob = p)
      try(fit.train <- fit.re(data[-test,]))
      fitted <- exists("fit.train")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    mc = lapply(fit.train, function(fit) {sum((as.numeric(predict(fit, newdata = data[test,])) - data[test,'sl'])^2)/g})
    adj = lapply(fit.train, function(fit) {sum((as.numeric(predict(fit, newdata = data)) - data[,'sl'])^2)/n})
    mc.ba = as.numeric(mc) + ws.mc - as.numeric(adj)
    list(mc = mc, mc.ba = mc.ba)
  }, mc.cores = 50)
  #return(CV.iter)
  models <- names(CV.iter[[1]]$mc)
  mci <- sapply(CV.iter, function(x) as.numeric(x$mc))#; return(CV.iter)
  mci.ba <- sapply(CV.iter, function(x) as.numeric(x$mc.ba))
  mc <- apply(mci,1,mean)#; return(MSE)
  mc.ba <- apply(mci.ba,1,mean)
  err <- apply(mci,1,sd)*sqrt(g)/sqrt(nrow(data))
  err.ba <- apply(mci.ba,1,sd)*sqrt(g)/sqrt(nrow(data))
  fsel <- apply(apply(mci,2, function(b) b == min(b)),1,sum)
  fsel.ba <- apply(apply(mci.ba,2, function(b) b == min(b)),1,sum)
  tibble(Model = models, 
         mc = mc, mc.ba = mc.ba, 
         err = err, err.ba = err.ba,
         Score = n*log(mc), Score.ba = n*log(mc.ba), 
         Delta.Score = Score - min(Score), Delta.Score.ba = Score.ba - min(Score.ba),
         fsel = fsel/sum(fsel), fsel.ba = fsel.ba/sum(fsel.ba))
}




# calculates AIC and BIC for conditional log-likelihood
IT.fit <- function(data){
  models <- fit.re(data)
  ll <- ll.cond(models)
  list(aic = -2*ll + 2*rep(c(0,1,1,2),3), bic = -2*ll + log(n)*rep(c(0,1,1,2),3))
}

# computes conditional log-likelihood for a list of nlme models
ll.cond <- function(fits) sapply(fits, function(m) -(n/2)*log(mean(as.numeric(m$residuals[,'haul'])^2)))

# generates a single bootstrap sample and returns IC
IT.boot.fit <- function(data){
  fitted <- F
  while (!fitted){
    data.boot <- sample_n(data,nrow(data), replace = T)
    models <- try(fit.re(data.boot))
    fitted <- !as.logical(sum(sapply(models,is.character)))
  }
  ll <- ll.cond(models)
  list(aic = -2*ll + 2*rep(c(0,1,1,2),3), bic = -2*ll + log(n)*rep(c(0,1,1,2),3))
  #list(aic = sapply(models, AIC), bic = sapply(models, BIC))
}

# frequency selection probabilities for AIC and BIC
IT.boot <- function(data, nboot){
  raw.boot <- pbmclapply(1:nboot, function(i) IT.boot.fit(data),mc.cores = 50)
  fsel.aic <- apply(sapply(raw.boot, function(b) b$aic == min(b$aic)),1,sum)
  fsel.bic <- apply(sapply(raw.boot, function(b) b$bic == min(b$bic)),1,sum)
  list(aic = fsel.aic/sum(fsel.aic), bic = fsel.bic/sum(fsel.bic))
}


}


#-------------
# Random Ideas
#-------------

if (F){
  

# random g CV
fe.rpcv <- function(iter, data, maxIter = 1000){
  n <- nrow(data)
  # p gives sampling probabilities proportionate to haul size (guarantees at least n.min random training samples in each haul)
  n.min <- 3
  p <- data %>% group_by(haul) %>% mutate(p = as.numeric(!1:n() %in% sample(n(),n.min))) %>% pull(p) 
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; it <- 0
    g <- sample(1:ceiling(n*(1 - (log(n) - 1)^(-1))),1)
    while(!fitted){
      it <- it + 1
      test <- sample(1:n, g, prob = p, replace = F)
      try(fit.train <- fit.re.fe(data[-test,]))
      fitted <- exists("fit.train")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/g
    fe.adj2 <- (sum(as.numeric(predict(fit.train[[1]], newdata = data)) - data[,'sl'])^2)/nrow(data)
    list(mse = rep(fe,g), adj2 = rep(as.numeric(fe.adj2),g))
    #re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,], level = 1)) - data[test,'sl'])^2)/g
    #c(fe= fe,re = re)
  }, mc.cores = 50)
  return(CV.iter)
}

# wildly variable, especially adj2 values

#CV.iter <- fe.rpcv(100000,fishData)
#MSE <- sapply(CV.iter, function(x) as.numeric(x$mse)) %>% simplify() %>% mean; MSE
#adj2 <- sapply(CV.iter, function(x) as.numeric(x$adj2)) %>% simplify() %>% mean; adj2

#ws.MSE <- fit.re.fe(fishData)$fe$m$resid()^2 %>% mean; ws.MSE
#BC <- ws.MSE - adj2; BC
#(n*log(ws.MSE) - n*log(MSE))/-2





# loocv per group (haul)
new.loopg <- function(iter, data, maxIter = 1000){
  n.group <- length(unique(data$haul))
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; iter <- 0
    while(!fitted){
      iter <- iter + 1
      test <- data %>% group_by(haul) %>% mutate(p = (1:n() %in% sample(n(),1))) %>% pull(p)
      try(fit.train <- fit.re.fe(data[-test,]))
      fitted <- exists("fit.train")
      if (iter > maxIter) stop("MaxIter exceeded")
    }
    fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/n.group
    re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,], level = 1)) - data[test,'sl'])^2)/n.group
    c(fe= fe,re = re)
  }, mc.cores = 50)
  #return(CV.iter)
  CV.iter <- sapply(CV.iter, function(x) x)#; return(CV.iter)
  MSE <- apply(CV.iter,1,mean)
  MSE.se <- apply(CV.iter,1,sd)*sqrt(n.group)/sqrt(nrow(data))
  fsel <- apply(apply(CV.iter,2, function(b) b == min(b)),1,sum)
  tibble(Model = rownames(CV.iter), MSE = MSE, MSE.se = MSE.se, Score = nrow(data)*log(MSE), 
         Delta.Score = Score - min(Score), fsel = fsel/sum(fsel))
}

#a <- new.loopg(500,fishData)
#a




# bias-adjusted loopg CV
adj.loopg <- function(iter, data, maxIter = 1000){
  n.group <- length(unique(data$haul))
  CV.iter <- pbmclapply(1:iter, function(i){
    fitted <- F; it <- 0
    while(!fitted){
      it <- it + 1
      test <- data %>% group_by(haul) %>% mutate(p = (1:n() %in% sample(n(),1))) %>% pull(p)
      try(fit.train <- fit.re.fe(data[-test,]))
      fitted <- exists("fit.train")
      if (it > maxIter) stop("MaxIter exceeded")
    }
    mse.fe <- sum((as.numeric(predict(fit.train[[1]], newdata = data[test,])) - data[test,'sl'])^2)/n.group
    adj2.fe <- (sum((as.numeric(predict(fit.train[[1]], newdata = data)) - data[,'sl'])^2)/nrow(data))/iter
    mse.re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,])) - data[test,'sl'])^2)/n.group
    adj2.re <- (sum((as.numeric(predict(fit.train[[2]], newdata = data)) - data[,'sl'])^2)/nrow(data))/iter
    #re <- sum((as.numeric(predict(fit.train[[2]], newdata = data[test,], level = 1)) - data[test,'sl'])^2)/g
    list(fe = c(mse.fe = mse.fe, adj2.fe = adj2.fe), re = c(mse.re= mse.re, adj2.re = adj2.re))
  }, mc.cores = 50)
  CV.iter <- sapply(CV.iter, function(x) x); return(CV.iter)
}


randomPartition <- function(n){
  smax <- ceiling(n*(1 - (log(n) - 1)^(-1)))
  n.this <- n
  i <- 1
  part <- c()
  while (n.this>0){
    s <- sample(1:min(smax,n.this),1)
    n.this <- n.this - s
    part <- c(part,s)
  }
  part
} 


} # end random ideas


#-------------
# OLD Stage 1 funs
#-------------

if(F){

  
  # no effects (ignore haul and sex), von B(vb), Gomperts(G) and logistic(lg)
  fit.ne <- function(data, sv = list(Linf = 210, K = 1, t0 = 0)){
    mods <- list()
    mods$vB = sl ~ Linf * (1 - exp(-(K * (age - t0))))
    mods$G = sl~ Linf * exp(-exp(-K * (age - t0)))
    mods$log = sl ~ Linf/(1 + exp(-K * (age - t0)))
    control = list(maxiter=10000,minFactor=1e-5,tol=1e-5)
    lapply(mods, function(m) nls(m, data = data, start = sv, control = control))
  }
  
  
  # fixed effects, haul and sex models, vb
  fit.vb.fe <- function(data, sv = c(Linf = 210, K = 0.3, t0 = -1)){
    mods <- c('0' = "0", K = "K",t = "t", Kt = "Kt")
    vb.fe <-  sv.vb <- list()
    vb.fe$'0' <- sl~Linf[haul]*(1-exp(-K*(age-t0)))
    vb.fe$K <- sl~Linf[haul]*(1-exp(-K[sex]*(age-t0)))
    vb.fe$t <- sl~Linf[haul]*(1-exp(-K*(age-t0[sex])))
    vb.fe$Kt <- sl~Linf[haul]*(1-exp(-K[sex]*(age-t0[sex])))
    
    sv.vb$'0' <- mapply(rep,sv,c(n.haul,1,1))
    sv.vb$K <- mapply(rep,sv,c(n.haul,2,1))
    sv.vb$t <- mapply(rep,sv,c(n.haul,1,2))
    sv.vb$Kt <- mapply(rep,sv,c(n.haul,2,2))
    
    lapply(mods, function(mod) nls(vb.fe[[mod]], data=data, start=sv.vb[[mod]]))
  }
  
  # fixed effects, haul and sex models, G
  fit.G.fe <- function(data, sv = c(Linf = 210, K = 0.3, t0 = -1)){
    mods <- c('0' = "0", K = "K",t = "t", Kt = "Kt")
    G.fe <-  sv.G <- list()
    G.fe$'0' <- sl~ Linf[haul] * exp(-exp(-K * (age - t0)))
    G.fe$'K' <- sl~ Linf[haul] * exp(-exp(-K[sex] * (age - t0)))
    G.fe$'t' <- sl~ Linf[haul] * exp(-exp(-K * (age - t0[sex])))
    G.fe$'Kt' <- sl~ Linf[haul] * exp(-exp(-K[sex] * (age - t0[sex])))
    
    sv.G$'0' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,1,1))
    sv.G$'K' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,2,1))
    sv.G$'t' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,1,2))
    sv.G$'Kt' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,2,2))
    
    lapply(mods, function(mod) nls(G.fe[[mod]], data=data, start=sv.G[[mod]]))
  }
  
  # fixed effects, haul and sex models, lg
  fit.lg.fe <- function(data){
    mods <- c('0' = "0", K = "K",t = "t", Kt = "Kt")
    lg.fe <-  sv.lg <- list()
    lg.fe$'0' <- sl ~ Linf[haul]/(1 + exp(-K * (age - t0)))
    lg.fe$'K' <- sl ~ Linf[haul]/(1 + exp(-K[sex] * (age - t0)))
    lg.fe$'t' <- sl ~ Linf[haul]/(1 + exp(-K * (age - t0[sex])))
    lg.fe$'Kt' <- sl ~ Linf[haul]/(1 + exp(-K[sex] * (age - t0[sex])))
    
    sv.lg$'0' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,1,1))
    sv.lg$'K' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,2,1))
    sv.lg$'t' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,1,2))
    sv.lg$'Kt' <- mapply(rep,c(Linf = 210, K = 0.3, t0 = -1),c(n.haul,2,2))
    
    lapply(mods, function(mod) nls(lg.fe[[mod]], data=data, start=sv.lg[[mod]]))
  }
  
  # fits no- and fixed- effects models for all three growth functions
  fit.ne.fe <- function(data) c(ne = fit.ne(data), vB = fit.vb.fe(data), G = fit.G.fe(data), lg = fit.lg.fe(data))
  
  
} # end old stage 1 funs


