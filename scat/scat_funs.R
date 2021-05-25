
prep_data <- function(data){
  data %>% as_tibble %>% 
    select(y = Species, Number, Length, Diameter, Taper, TI, Mass, CN, Location, ropey, segmented) %>% 
    rename_with(str_to_lower) %>% 
    mutate(across(c(mass,cn),log)) %>% # following original data publication
    mutate(number = as.double(number)) %>% 
    drop_na() %>% 
    mutate(across(where(is_double), ~ .x %>% scale %>% as.numeric)) %>% # centre and scale predictors
    mutate(y = fct_collapse(y, canid = c("coyote","gray_fox"), felid = c("bobcat")) %>% 
             fct_relevel("canid","felid"),
           location = fct_collapse(location, mid = c("middle"), edge = c("edge","off_edge")) %>% 
             factor(labels = c(0,1))) %>% 
    mutate(across(c(ropey,segmented), factor))
}


step_glm <- function(vars, cv_data, metric = "tss"){
  if(metric == "log_density") fit_fun <- fit_ll_glm
  if(metric == "tss") fit_fun <- fit_confusion_glm
  
  sel <- c()
  metric_step <- list()
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
  
  for(step in 1:length(vars)){
    message(paste("Step:",step))
    newfits = pbmclapply(setdiff(vars,sel), function(var){
      fit_fun(cv_data, make_formula(c(sel,var)))
    }, mc.cores = length(setdiff(vars,sel)))
    best = newfits %>% map("metric") %>% map_dbl(mean) %>% which.max()
    sel[step] <- setdiff(vars,sel)[best]
    message(paste("Selected:",sel[step]))
    metric_step[[step]] <- newfits[[best]]
  }
  list(sel = sel, metric_step = metric_step)
}


fit_confusion_glm <- function(cv_data, formula = NULL){
  get_confusion_matrix <- function(form, split, thresh){
    glm(form, data = analysis(split), family = binomial()) %>% 
      {tibble(pred = predict(.,newdata = assessment(split), type = "response") >= thresh, 
              y = assessment(split)$y != levels(assessment(split)$y)[1])} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  if(!is.null(formula)) form <- list(formula) else form <- cv_data$form_glm
  with(cv_data, mapply(get_confusion_matrix, form, splits, thresh_glm, SIMPLIFY = T)) %>% t %>% 
    as_tibble %>% 
    bind_cols(rep = cv_data$id, .) %>% 
    group_by(rep) %>% 
    summarise(metric = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% tss) %>% 
    mutate(rep = 1:n())
}

fit_ll_glm <- function(cv_data, formula){
  calc_pred_ll <- function(form, split){
    fit <- glm(form, data = analysis(split), family = binomial())
    tibble(y_pred = predict(fit, newdata = assessment(split),type = "response"),
           y_data = as.numeric(assessment(split)$y != levels(assessment(split)$y)[1]),
           log_p = log(abs(y_data - 1 + y_pred))) %>% pull(log_p) %>% sum
  }
  with(cv_data, mapply(calc_pred_ll, list(formula), splits, SIMPLIFY = T)) %>% 
    bind_cols(rep = cv_data$id, ll = .) %>% 
    group_by(rep) %>% 
    summarise(metric = sum(ll)) %>% 
    mutate(rep = 1:n())
}


tss <- function (cm) {
  sens <- cm[1,1]/(cm[2,1] + cm[1,1])
  spec <- cm[2,2]/(cm[1,2] + cm[2,2])
  sens + spec - 1
}


tss_thresh <- function(tbl){
  probs<- tbl$probs
  y <- tbl$y
  probs_to_tss <- function(thresh){
    pred  = as.numeric(probs >= thresh) %>% as.numeric %>% factor(levels = c(1,0))
    y = y %>% factor(levels = c(1,0))
    tss(table(pred,y))
  }
  threshold = (1:100)/100
  tss =  threshold %>% map_dbl(probs_to_tss)
  tibble(threshold,tss) %>% arrange(-tss) %>% slice(1) %>% as_vector
}


# cv_data: nested cv object from rsamples
tune_rf <- function(cv_data, mtry_vec, ntree){
  pbmclapply(cv_data$inner_resamples, function(sample){
    lapply(mtry_vec, function(m){
      sample[["splits"]] %>% map(
        ~ randomForest(y ~ ., mtry = m, data = analysis(.x), ntree = ntree) %>% 
          predict(newdata = assessment(.x), type = "prob") %>% 
          {tibble(probs = .[,"felid"], y = as.numeric(assessment(.x)$y =="felid"))}
      ) %>% bind_rows %>% tss_thresh
    }) %>% bind_rows %>% mutate(mtry = mtry_vec) %>% arrange(-tss) %>% slice(1)
  }, mc.cores = MAX_CORES) %>% bind_rows
}


tune_glm <- function(cv_data){
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
    pbmclapply(cv_data_2$inner_resamples, function(sample){
    sel <- thresh <- c()
    metric_step = list()
    
    for(step in 1:length(vars)){
      message(paste("Step:",step))
      newfits = lapply(setdiff(vars,sel), function(var){
        sample[["splits"]] %>% map(
          ~ glm(make_formula(c(sel,var)), data = analysis(.x), family = binomial()) %>% 
            {tibble(probs = predict(.,newdata = assessment(.x), type = "response"), 
                    y = as.numeric(assessment(.x)$y != levels(assessment(.x)$y)[1]))}
        ) %>% bind_rows %>% tss_thresh
      })
      best = newfits %>% map("tss") %>% map_dbl(mean) %>% which.max()
      thresh[step] = newfits[[best]][["threshold"]]
      sel[step] = setdiff(vars,sel)[best]
      message(paste("Selected:",sel[step]))
      metric_step[[step]] = newfits[[best]]
    }
    best_model = metric_step %>% map_dbl("tss") %>% which.max
    list(form_glm = make_formula(sel[1:best_model]), 
         thresh_glm = thresh[best_model])
  }, mc.cores = MAX_CORES)
}


fit_rf <- function(cv_data, ntree = 500, type = c("best","all")){
  make_formula = function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + ")) %>% as.formula()
  type = type[1]
  get_confusion_matrix <- function(split, thresh, mtry, ntree){
    success = levels(assessment(split)$y)[2]
    if(type == "best"){ # select the most important vars up to mtry terms
      fit_init = randomForest(y ~ ., mtry = mtry, data = analysis(split), ntree = ntree)
      form = fit_init %>% importance() %>% as.data.frame %>% 
        arrange(-MeanDecreaseGini) %>% slice(1:mtry) %>% row.names() %>% make_formula
    } else form = formula(y~.)
    randomForest(form, mtry = mtry, data = analysis(split), ntree = ntree) %>% 
      {tibble(pred = predict(., newdata = assessment(split), type = "prob")[,success] >= thresh,
              y = assessment(split)$y == success)} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  with(cv_data, pbmcmapply(get_confusion_matrix, splits, thresh_rf, mtry_rf, list(ntree), mc.cores = MAX_CORES)) %>% 
    simplify %>% t %>% 
    as_tibble %>% 
    mutate(rep = 1:n()) %>% 
    group_by(rep) %>% 
    summarise(metric = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% tss)
}

make_plot_data <- function(metric_data, levels){
  best_model <- metric_data %>% map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% map_dbl(mean),
         se = metric_data %>% map_dbl(~ sd(.x)/sqrt(length(.x))),
         se_diff = metric_data %>% map(~ .x - metric_data[[best_model]]) %>% map_dbl(~ sd(.x)/sqrt(length(.x))),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}


plot_model_comparisons <- function(plot_data, se_type = c("se_mod", "se_diff", "se")){
  plot_data <- plot_data %>% mutate(se = plot_data[[se_type[1]]])
  plot_data %>% 
    ggplot(aes(model)) +
    geom_point(aes(y = metric)) +
    geom_linerange(aes(ymin = metric - se, ymax = metric + se)) +
    theme_bw() 
}




# OLD STUFF

if(F){






step_glm <- function(vars, cv_data, metric = "log_density"){
  if(metric == "log_density") fit_fun <- fit_ll_glm
  if(metric == "tss") fit_fun <- fit_confusion_glm
  
  sel <- c()
  metric_step <- list()
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
  
  for(step in 1:length(vars)){
    message(paste("Step:",step))
    newfits = pbmclapply(setdiff(vars,sel), function(var){
      fit_fun(make_formula(c(sel,var)),cv_data)
    }, mc.cores = length(setdiff(vars,sel)))
    best = newfits %>% map("metric") %>% map_dbl(mean) %>% which.max()
    sel[step] <- setdiff(vars,sel)[best]
    message(paste("Selected:",sel[step]))
    metric_step[[step]] <- newfits[[best]]
  }
  list(sel = sel, metric_step = metric_step)
}




step_glm_2 <- function(vars, cv_data, metric = "log_density"){
  if(metric == "log_density") fit_fun <- fit_ll_glm
  if(metric == "tss") fit_fun <- fit_confusion_glm
  
  sel <- c()
  metric_step <- list()
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
  
  for(step in 1:length(vars)){
    message(paste("Step:",step))
    newfits = lapply(setdiff(vars,sel), function(var){
      fit_fun(make_formula(c(sel,var)),cv_data)
    })
    best = newfits %>% map("metric") %>% map_dbl(mean) %>% which.max()
    sel[step] <- setdiff(vars,sel)[best]
    message(paste("Selected:",sel[step]))
    metric_step[[step]] <- newfits[[best]]
  }
  list(sel = sel, metric_step = metric_step %>% map_dbl("metric"))
}







tune_tss_mtry_old <- function(data, mtry_vec, K = 10, ntree = 500){
  folds = sample(rep(1:K, length.out = nrow(data)))
  sapply(mtry_vec, function(m){
    1:K %>% 
      map(~ randomForest(y ~ ., mtry = m, data = data[folds!=.x,], ntree = ntree) %>% 
            predict(newdata = data[folds==.x,], type = "prob") %>% 
            {tibble(probs = .[,"felid"], y = as.numeric((data[folds==.x,"y"]) =="felid"))}) %>% 
      bind_rows %>% tss_thresh 
  }) %>% t %>% as.data.frame %>% as_tibble %>% mutate(mtry = mtry_vec)
}

}
           