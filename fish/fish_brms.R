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
library(fishmethods)
library(FlexParamCurve)
library(RColorBrewer)
library(brms)

rm(list=ls())

options(brms.file_refit = "on_change")

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
min.size <- 10
fishData <- fishData %>% group_by(haul) %>% filter(n()>= min.size) %>% ungroup
fishData <- fishData %>% filter(haul != "TBD930066") # remove outlier
fishData$haul <-  factor(fishData$haul)

# MLE fits for initial values and priors
f <- fit.re(fishData)
f$log.Kt %>% coef %>% map_dbl(mean)


run_date <- "2021_10_01"
sv_dir <- paste0("fits_",run_date)

#----------
# vB models
#----------

f.vB.0 <-  bf(sl ~ (150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
f.vB.K <-  bf(sl ~ (150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ sex, t0 ~ 1, Linf ~ 1 + 1|haul, nl = TRUE)
f.vB.t <-  bf(sl ~ (150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ sex, Linf ~ 1 + 1|haul, nl = TRUE)
f.vB.Kt <-  bf(sl ~ (150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ sex, t0 ~ sex, Linf ~ 1 + 1|haul, nl = TRUE)
f.vB.L <-  bf(sl ~ (150 + Linf) * (1 - exp(-1*((1 + K) * (age - t0)))), K ~ 1, t0 ~ 1, Linf ~ sex + (1|haul), nl = TRUE)
# N.B. vB.L does not improve over vB.0 or vB.t

pr.vB <- 
  prior(normal(0, 1), nlpar = "K") +  
  prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
  prior(normal(0, 1), nlpar = "t0") +
  prior(normal(0, 3), nlpar = "Linf")

pr.vB.K <- pr.vB + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") 
pr.vB.t <- pr.vB + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
pr.vB.Kt <- pr.vB.K + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
pr.vB.L <- pr.vB + prior(normal(0, 0.2), nlpar = "Linf", coef = "sex2") 

m.vB.0 <- brm(f.vB.0, data = fishData, prior = pr.vB, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.0"))
m.vB.K <- brm(f.vB.K, data = fishData, prior = pr.vB.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.K"))
m.vB.t <- brm(f.vB.t, data = fishData, prior = pr.vB.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.t"))
m.vB.Kt <- brm(f.vB.Kt, data = fishData, prior = pr.vB.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.Kt"))
m.vB.L <- brm(f.vB.L, data = fishData, prior = pr.vB.L, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.vB.L"))


m.vB <- loo::nlist(m.vB.0, m.vB.K, m.vB.t, m.vB.Kt)
loo.vB <- m.vB %>% map(loo, cores = 10)
loo.vB %>% loo_compare

#------------

#-----------------
# Gomperts models
#-----------------

f.G.0 <-  bf(sl~ (151 + Linf) * exp(-exp(-(1.25 + K) * (age - t0))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
f.G.K <-  bf(sl~ (151 + Linf) * exp(-exp(-(1.25 + K) * (age - t0))), K ~ sex, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
f.G.t <-  bf(sl~ (151 + Linf) * exp(-exp(-(1.25 + K) * (age - t0))), K ~ 1, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
f.G.Kt <-  bf(sl~ (151 + Linf) * exp(-exp(-(1.25 + K) * (age - t0))), K ~ sex, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)

pr.G <- 
  prior(normal(0, 1), nlpar = "K") +  
  prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
  prior(normal(0, 1), nlpar = "t0") +
  prior(normal(0, 3), nlpar = "Linf")

pr.G.K <- pr.G + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") 
pr.G.t <- pr.G + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
pr.G.Kt <- pr.G.K + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 

m.G.0 <- brm(f.G.0, data = fishData, prior = pr.G, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.0"))
m.G.K <- brm(f.G.K, data = fishData, prior = pr.G.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.K"))
m.G.t <- brm(f.G.t, data = fishData, prior = pr.G.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.t"))
m.G.Kt <- brm(f.G.Kt, data = fishData, prior = pr.G.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.G.Kt"))

m.G <- loo::nlist(m.G.0, m.G.K, m.G.t, m.G.Kt)
loo.G <- m.G %>% map(loo, cores = 20)
loo.G %>% loo_compare

#------------


#-----------------
# logistic models
#-----------------

f.log.0 <-  bf(sl ~ (151 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
f.log.K <-  bf(sl ~ (151 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ 1, Linf ~ 1 + (1|haul), nl = TRUE)
f.log.t <-  bf(sl ~ (151 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ 1, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)
f.log.Kt <-  bf(sl ~ (151 + Linf)/(1 + exp(-(1.4 + K) * (age - t0))), K ~ sex, t0 ~ sex, Linf ~ 1 + (1|haul), nl = TRUE)

pr.log <- 
  prior(normal(0, 1), nlpar = "K") +  
  prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
  prior(normal(0, 1), nlpar = "t0") +
  prior(normal(0, 3), nlpar = "Linf")

pr.log.K <- pr.log + prior(normal(0, 0.2), nlpar = "K", coef = "sex2") 
pr.log.t <- pr.log + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 
pr.log.Kt <- pr.log.K + prior(normal(0, 0.2), nlpar = "t0", coef = "sex2") 

m.log.0 <- brm(f.log.0, data = fishData, prior = pr.log, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.0"))
m.log.K <- brm(f.log.K, data = fishData, prior = pr.log.K, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.K"))
m.log.t <- brm(f.log.t, data = fishData, prior = pr.log.t, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.t"))
m.log.Kt <- brm(f.log.Kt, data = fishData, prior = pr.log.Kt, cores = 4, iter = 4000, file = paste0(sv_dir,"/m.log.Kt"))

m.log <- loo::nlist(m.log.0, m.log.K, m.log.t, m.log.Kt)
loo.log <- m.log %>% map(loo, cores = 10)
loo.log %>% loo_compare

#------------

loo.m <- c(loo.vB,loo.G,loo.log)
#saveRDS(loo.m, "fits_2021_10_01/loo.m.rds")

make_plot_data <- function(metric_data, levels = names(metric_data)){
  best_model <- metric_data %>% map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% map_dbl(mean),
         metric_diff = metric - max(metric),
         se = metric_data %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_diff = metric_data %>% map(~ .x - metric_data[[best_model]]) %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}

loo_df <- make_plot_data(loo_pointwise)

loo_df <- loo_df %>% 
  as_tibble() %>% 
  mutate(across(everything(), as.numeric)) %>% 
  mutate(model = names(loo.m)) %>% 
  mutate(fun = c(rep("vB",4),rep("G",4),rep("log",4)), 
         fe = rep(c("0","K","t","Kt"),3)) %>% 
  mutate(fe = factor(fe, levels = c("0","K","t","Kt")))
  
fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")

loo_df

# faceted plot
loo_df %>% 
  ggplot(aes(x = fe)) +
  #geom_line(aes(y = metric_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = metric_diff, col = fun), size = 2, show.legend = F) +
  geom_linerange(aes(ymin = metric_diff - se_mod, ymax = metric_diff + se_mod, col = fun), show.legend = F) +
  theme_classic() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90")) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = metric_diff), shape = 1, size = 6, col = "black", data = ~ .x %>% filter(model == "m.vB.0")) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  labs(x = NULL, subtitle = "Pinfish PSIS-LOO estimates", y = expression(Delta*"ELPD"))

ggsave("plots/pinfish_modsel_2021_10_01.pdf", width = 90, height = 80, units = "mm", device = cairo_pdf()); dev.off()  

