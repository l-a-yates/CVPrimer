#----------------------------------------------------------------------------
#
# Cross validation for model selection: a primer with examples from ecology.
#
# Code for example 2: Pinfish growth models
#
# Authors: L.A.Yates, ...
#
#----------------------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(fishmethods)
library(FlexParamCurve)
library(RColorBrewer)

source("fish_functions_flex.R")
source("fish_functions_nlme.R")

MAX_CORES <- 30


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

# constants
#n.haul <- length(unique(fishData$haul))
#n <- nrow(fishData)

#-----------------------------------------------------
# fit 4 growth functions for population-level analysis
#-----------------------------------------------------
if(F){
  # 4 par model: Richards; 3 par models: von Bertanlanffy, logistic, Gompertz
  m <- fitGrowthAll(fishData)
  # N.B. it is hard to fit the Richards model, nls only converges when using a large tolerance
  m %>% map_dbl(c("convInfo","finTol")) 
  # compare within-sample sum of least squares 
  # N.B. we can reject Richards immediately since it is the most complex and the worse fitting model
  m %>% map(~ .x %>% resid %>% {.^2} %>% sum) %>% unlist %>% {.-min(.)}
}

#---------------------------------------------------
# fit 3 growth functions as mixed effects models
#---------------------------------------------------

# a quick within sample fit
f1 <- fit.re(fishData)

# marginal AIC only (random effects are (erroneously) ignored in parameter count!)
# i.e. numbers of parameters = number of fixed effects; check: YES
f1 %>% map_dbl(AIC) %>% {. - min(.)} 
p_fe <- f1 %>% {map_dbl(.,AIC)*0.5 + map_dbl(.,logLik)}; p_fe

# cross validation estimates
f_cv_loo <- loocv(fishData, fit.re)  # 15 seconds wth 30 cores
f_cv_10 <- kfcv(fishData, fit.re, 10, bias.adjust = T) # 2 seconds with 10 cores
f_cv_10_nba <- kfcv(fishData, fit.re, 10, bias.adjust = F) # 2 seconds with 10 cores

# compare p_loo differences to the number of fixed effects
tibble(p_loo_dif = f_cv_loo$p_loo %>% {. - min(.)}, p_fe_diff = p_fe - 5)

# compare bias-adjusted loo estimates (almost identical)
f_cv_loo_nba <- loocv(fishData, fit.re, bias.adjust = F) # 15 seconds wth 30 cores
(f_cv_loo[,-1] - f_cv_loo_nba[,-1]) %>% bind_cols(model = f_cv_loo$model, .)


# prepares returned CV fits for plotting
prep_df <- function(df){
  df %>% 
    mutate(fun = c(rep("vB",4),rep("G",4),rep("log",4)), 
           fe = rep(c("0","K","t","Kt"),3)) %>% 
    mutate(fe = factor(fe, levels = c("0","K","t","Kt")))
}

fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")

# colors
#display.brewer.all(colorblindFriendly = TRUE)
#display.brewer.pal(n = 8, name = 'Dark2')

# non-faceted plot
f_cv_loo %>% 
  prep_df %>% 
  mutate(model = str_replace(model,fixed("."),"|")) %>% 
  ggplot(aes(x = model)) +
  geom_line(aes(y = eld_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = eld_diff, col = fun), size = 3, show.legend = F) +
  geom_linerange(aes(ymin = eld_diff - se_mod, ymax = eld_diff + se_mod, col = fun), show.legend = F) +
  theme_bw() +
  scale_color_brewer(type = "div", palette = "Dark2")+
  geom_point(aes(y = eld_diff), shape = 1, size = 6, col = "grey40", data = ~ .x %>% filter(model == "vB|0")) +
  labs(x = NULL, title = "LOO estimates for fish growth models", y = expression(Delta*"ELPD"))

# faceted plot
f_cv_loo %>% 
  prep_df %>% 
  #mutate(model = str_replace(model,fixed("."),"|")) %>% 
  ggplot(aes(x = fe)) +
  geom_line(aes(y = eld_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = eld_diff, col = fun), size = 3, show.legend = F) +
  geom_linerange(aes(ymin = eld_diff - se_mod, ymax = eld_diff + se_mod, col = fun), show.legend = F) +
  theme_bw() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 11), strip.background.x = element_rect(linetype = 0, fill = "grey90")) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = eld_diff), shape = 1, size = 6, col = "grey40", data = ~ .x %>% filter(model == "vB.0")) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  labs(x = NULL, title = "LOO estimates for fish-growth models", y = expression(Delta*"ELPD"))

#ggsave("plots/pinfish_modsel_loo_v2.pdf", width = 180, height = 100, units = "mm", device = cairo_pdf()); dev.off()

# faceted plot comparing to score estimators
alt_data <- f_cv_10 %>% mutate(score = eld_diff, se = se_mod)
alt_data2 <- f_cv_10 %>% mutate(score = eld_diff, se = se_mod)

# faceted plot
f_cv_loo %>% 
  prep_df %>% 
  mutate(score = eld_diff, se = se_mod) %>% 
  ggplot(aes(x = fe)) +
  geom_line(aes(y = score, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_line(aes(y = score, group = fun), col = "grey40", lty = "dashed", show.legend = F,
            alpha = 0.1, position = position_nudge(0.15), data = alt_data %>% prep_df) +
  #geom_line(aes(y = score, group = fun), col = "grey40", lty = "dashed", show.legend = F,
  #          alpha = 0.4, position = position_nudge(0.15), data = alt_data2 %>% prep_df) +
  
  geom_point(aes(y = score, col = fun), size = 3, show.legend = F) +
  geom_point(aes(y = score, col = fun), alpha = 0.1, size = 3, show.legend = F, 
             position = position_nudge(0.15), data = alt_data %>% prep_df) +
  geom_point(aes(y = score, col = fun), alpha = 0.1, size = 3, show.legend = F, 
             position = position_nudge(0.15), data = alt_data2 %>% prep_df) +
  
  geom_linerange(aes(ymin = score - se, ymax = score + se, col = fun), show.legend = F,
                 alpha = 0.1, position = position_nudge(0.15), data = alt_data %>% prep_df) +
  geom_linerange(aes(ymin = score - se, ymax = score + se, col = fun), show.legend = F,
                 alpha = 0.1, position = position_nudge(0.15), data = alt_data2 %>% prep_df) +
  geom_linerange(aes(ymin = score - se, ymax = score + se, col = fun), show.legend = F) +
  theme_bw() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank() ) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = score), shape = 1, size = 6, col = "grey40", data = ~ .x %>% filter(model == "vB.0")) +
  labs(x = NULL, title = "CV estimates for fish growth models", y = "score",
       subtitle = "A comparison between LOO and bias-adjusted 10-fold CV")

#ggsave("pinfish_modsel_compare_ba_v1.pdf", width = 210, height = 160, units = "mm", device = cairo_pdf()); dev.off()


#----------------------------------------------------------
# Compare best model using haul as a random vs fixed effect
#----------------------------------------------------------

f2_loo <- loocv(fishData, fit.re.fe.0)

f2_loo %>% select(model, eld_diff, se_mod, p_loo) %>% 
  mutate(across(where(is.numeric), round, digits = 3)) %>% 
  kable(format = "latex", booktabs = TRUE)

f2_loo %>% 
  mutate(model = factor(model, levels = c("re","fe"), labels = c("Random","Fixed"))) %>% 
  ggplot(aes(model)) +
  geom_point(aes(y = eld_diff)) + 
  geom_linerange(aes(ymin = eld_diff - se_diff, ymax = eld_diff + se_diff)) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.title = element_blank()) 


display.brewer.pal(n = 8, name = 'Dark2')

f2_ws <- fit.re.fe.0(fishData)
f2_ws %>% map_dbl(AIC) %>% {. - min(.)} # delta AIC: completely wrong!

eff <- tibble(Fixed = f2_ws$fe %>% coef %>% {.[str_starts(names(.),"Linf")]},
       Random = f2_ws$re %>% coef %>% pull(Linf)) %>% 
  arrange(Random) %>% 
  mutate(haul = factor(1:length(Fixed)))

eff %>% 
  pivot_longer(1:2, names_to = "type", values_to = "L") %>% 
  ggplot(aes(x = haul)) +
  geom_point(aes(y = L, col = type)) +
  geom_hline(aes(yintercept = mean(eff$Random)), col = "grey70", lty = "dashed")+
  #ylim(c(NA,NA)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.title = element_blank()) +
  scale_color_brewer(type = "div", palette = "Dark2") +
  labs(title = expression("Haul-level estimates of"~L[0]), 
       subtitle = "Shrinkage of random with respect to fixed estimates",
       y = expression(L[0]), x = "Haul")


#ggsave("pinfish_shrinkage_plot_v1.pdf", width = 180, height = 110, units = "mm", device = cairo_pdf())


#----------
# Data plot
#----------

# colors
#display.brewer.all(colorblindFriendly = TRUE)
#display.brewer.pal(n = 8, name = 'Dark2')

age_lims <- fishData %>% {c(min = min(.$age), max = max(.$age))}
xs <- seq(0.5,6.3,0.02)

hauls_ordered <- f2_ws$re %>% coef %>% arrange(Linf) %>% rownames
Linfs <- f2_ws$re %>% coef %>% arrange(Linf) %>% pull(Linf)
names(Linfs) <- hauls_ordered

curve_data <- hauls_ordered %>% 
  map_dfr(~ tibble(haul = .x, 
                   age = xs,
                   l = predict(f2_ws$re, newdata = tibble(haul = .x,age = xs)),
                   Linf = Linfs[.x]))

point_data <- fishData %>% mutate(Linf = Linfs[as.character(haul)])

coefs_fixed <- f2_ws$re$coefficients$fixed
mean_model <- function(age) coefs_fixed["Linf"]* (1 - exp(-(coefs_fixed["K"]* (age - coefs_fixed["t0"]))))

curve_data %>% 
  ggplot(aes(age)) +
  geom_line(aes(y = l, group = haul, col = Linf), alpha = 0.4) +
  geom_line(aes(y = l), size = 1.5, lty = "longdash", data = tibble(age = xs, l = mean_model(age))) +
  geom_point(aes(y = sl, col = Linf), alpha = 0.7, data = point_data) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Age (years)", y ="Length (mm)", title = "Pinfish growth data")
  
#ggsave("pinfish_data_plot_v1.pdf", width = 180, height = 110, units = "mm", device = cairo_pdf())



stop("whoa!")

# Old stuff
if(F){
  
  # faceted model-selection plot
  strip_labels = c(`0` = "No fixed effects",K = "K ~ sex", t = "t ~ sex", Kt = "K + t ~ sex")
  f1_eld %>% 
    prep_df %>% 
    ggplot(aes(x = fun)) +
    geom_point(aes(y = eld_diff, col = fun), size =2) +
    geom_linerange(aes(ymin = eld_diff - se_mod, ymax = eld_diff + se_mod, col = fun)) +
    facet_wrap(~fe,nrow = 1, strip.position = "bottom", labeller = labeller(fe = strip_labels)) +
    theme_bw() +
    theme(strip.placement = "outside", panel.border = element_blank(), 
          panel.grid = element_blank(), axis.ticks.x = element_blank() ) +
    geom_point(aes(y = eld_diff), shape = 1, size = 5, col = "darkblue", data = ~ .x %>% filter(model == "vB.0")) +
    labs(x = NULL, title = "LOO estimates for fish growth models", y = "eld_diff")
  
  
  model_levels <- c("log.0","G.0","vB.0","log.K","G.K","vB.K","log.t","G.t","vB.t","log.Kt","G.Kt","vB.Kt")
  f1_cv %>% 
    mutate(model = factor(model, levels = model_levels)) %>% 
    mutate(loo = mse, se = se_mod) %>% 
    ggplot(aes(model)) +
    geom_point(aes(y = loo)) +
    geom_linerange(aes(ymin = loo -se, ymax = loo + se)) +
    labs(title = "LOO estimates for fish growth models") +
    geom_point(aes(y = loo), shape = 1, size = 4, col = "darkblue", data = ~ .x %>% filter(model == "vB.0")) +
    theme_classic()
  
  f1_cv %>% 
    arrange(p_loo) %>% 
    mutate(model = factor(model, levels = model)) %>% 
    mutate(loo = delta_mse, se = se_mod) %>% 
    ggplot(aes(model)) +
    geom_point(aes(y = loo)) +
    geom_linerange(aes(ymin = loo -se, ymax = loo + se)) +
    labs(title = "LOO estimates for fish-growth models",
         subtitle = "Models ordered by estimated complexity") +
    geom_point(aes(y = loo), shape = 1, size = 4, col = "darkblue", data = ~ .x %>% filter(model == "vB.0")) +
    theme_classic()
  
  f1_cv
  
  
  
}