
# using the FlexParamCurve package
fitModel <- function(type, data){
  if(!type %in% c("Gompertz","logistic","vonB","Richards")) stop("select valid type")
  modpar(fishData$age, fishData$sl, force4par = T, pn.options = "opts", verbose = F,suppress.text = T)
  if(type == "Richards"){
    notFit <- T
    tol = 20
    print("Fitting Richards model...")
    while(notFit){
      fit <-try(nls(sl ~ SSposnegRichards(age,Asym=Asym,K=K,Infl=Infl,M=M, modno = 12, pn.options = "opts"),
                    data = data,
                    control = nls.control(maxiter = 10, tol = tol*1e-3)), silent = T)
      if(!attr(fit, "class")=="try-error") notFit <- F else tol <- tol + 1
      if(tol > 50) stop("exceeded tolerance limit")
    }
    return(fit)
  }
  if(type == "Gompertz") M <- 0.1
  if(type == "logistic") M <- 1
  if(type == "vonB") M <- -0.3
  change.pnparameters(M=M, pn.options = "opts")
  print(paste("Fitting",type,"model..."))
  notFit <- T
  tol = 1
  while(notFit){
    fit <-try(nls(sl ~ SSposnegRichards(age,Asym=Asym,K=K,Infl=Infl, modno = 32, pn.options = "opts"),
                  data = data,
                  control = nls.control(maxiter = 50, tol = tol*1e-6)), silent = T)
    if(!attr(fit, "class")=="try-error") notFit <- F else tol <- tol + 1
    if(tol > 50) stop("exceeded tolerance limit")
  }
  return(fit)
}

# using the FlexParamCurve package
fitGrowthAll <- function(data){
  c(rich = "Richards", vonb = "vonB", logi = "logistic", gomp = "Gompertz") %>% 
    lapply(fitModel, data = data)
}



