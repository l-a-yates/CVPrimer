
prep_data <- function(data){
  data %>% as_tibble %>% 
    select(y = Species, Number, Length, Diameter, Taper, TI, Mass, CN, Location, ropey, segmented) %>% 
    rename_with(str_to_lower) %>% 
    mutate(across(c(mass,cn),log)) %>% # following original data publication
    drop_na() %>% 
    mutate(across(where(is_double), ~ .x %>% scale %>% as.numeric)) %>% # centre and scale predictors
    mutate(y = fct_collapse(y, canid = c("coyote","gray_fox"), felid = c("bobcat")) %>% 
             fct_relevel("canid","felid"),
           location = fct_collapse(location, mid = c("middle"), edge = c("edge","off_edge")) %>% 
             factor(labels = c(0,1)))
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


tune_tss_mtry <- function(data, mtry_vec, K = 10, ntree = 500){
  folds = sample(rep(1:K, length.out = nrow(data)))
  sapply(mtry_vec, function(m){
    1:K %>% 
      map(~ randomForest(y ~ ., mtry = m, data = data[folds!=.x,], ntree = ntree) %>% 
            predict(newdata = data[folds==.x,], type = "prob") %>% 
            {tibble(probs = .[,"felid"], y = as.numeric((data[folds==.x,"y"]) =="felid"))}) %>% 
      bind_rows %>% tss_thresh 
  }) %>% t %>% as.data.frame %>% as_tibble %>% mutate(mtry = mtry_vec)
}

fit_confusion_rf <- function(cv_data, ntree = 500){
  get_confusion_matrix <- function(split, thresh, mtry, ntree){
    success = levels(assessment(split)$y)[2]
    randomForest(y ~ ., mtry = mtry, data = analysis(split), ntree = ntree) %>% 
      {tibble(pred = predict(., newdata = assessment(split), type = "prob")[,success] >= thresh,
              y = assessment(split)$y == success)} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  with(cv_data, mapply(get_confusion_matrix, splits, thresh, mtry, list(ntree), SIMPLIFY = T)) %>% t %>% 
    as_tibble %>% 
    bind_cols(rep = cv_data$id, .) %>% 
    group_by(rep) %>% 
    summarise(tss = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% tss) %>% 
    mutate(rep = 1:n())
}


fit_confusion_glm <- function(formula, cv_data){
  get_confusion_matrix <- function(form, split, thresh){
    glm(form, data = analysis(split), family = binomial()) %>% 
      {tibble(pred = predict(.,newdata = assessment(split), type = "response") >= thresh, 
              y = assessment(split)$y != levels(assessment(split)$y)[1])} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  with(cv_data, mapply(get_confusion_matrix, list(formula), splits, thresh, SIMPLIFY = T)) %>% t %>% 
    as_tibble %>% 
    bind_cols(rep = cv_data$id, .) %>% 
    group_by(rep) %>% 
    summarise(tss = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% tss) %>% 
    mutate(rep = 1:n())
}
