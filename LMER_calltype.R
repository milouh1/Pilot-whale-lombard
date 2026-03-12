## Mixed-effects analysis: Call type
## =========================================================

## 0) Load packages 
library(lme4)
library(lmerTest)
library(DHARMa)
library(dplyr)
library(ggplot2)
library(car)
library(performance)
library(emmeans)
library(ggeffects)
library(sjPlot)

emm_options(lmer.df = "satterthwaite")  

## 1) Load data 
pwdata <- read.csv(
  "~/Desktop/stats_PW/Hegeman_Tarifa_PW_Parameters_2025_05_07.csv",
  stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

# Keep SNR >= 6
pwdata <- pwdata %>% filter(SNR >= 6)

# Exclude individuals (TAG) with < 20 calls
pwdata <- pwdata %>%
  group_by(TAG) %>%
  filter(n() >= 20) %>%
  ungroup()

# Factors
pwdata$CALL_TYPE <- factor(pwdata$CALL_TYPE)
pwdata$TAG       <- factor(pwdata$TAG)

#  "sk" as reference 
if ("sk" %in% levels(pwdata$CALL_TYPE)) {
  pwdata$CALL_TYPE <- relevel(pwdata$CALL_TYPE, ref = "sk")
}

# Within-tag centering of NL 
pwdata <- pwdata %>%
  group_by(TAG) %>%
  mutate(
    NL_bar = mean(NL, na.rm = TRUE),    
    NL_c   = NL - NL_bar                
  ) %>%
  ungroup()

## 2) Fitting the models
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 60000))

m0_w <- lmer(AOL_RMS ~ 1 + (1 + NL_c | TAG),
             data = pwdata, REML = FALSE, control = ctrl)

m1_w <- lmer(AOL_RMS ~ NL_c + (1 + NL_c | TAG),
             data = pwdata, REML = FALSE, control = ctrl)

m2_w <- lmer(AOL_RMS ~ NL_c + CALL_TYPE + (1 + NL_c | TAG),
             data = pwdata, REML = FALSE, control = ctrl)

## Full model
m3_w <- lmer(AOL_RMS ~ NL_c * CALL_TYPE + (1 + NL_c | TAG),
             data = pwdata, REML = FALSE, control = ctrl)   

# Likelihood-ratio tests
anova(m0_w, m1_w)
anova(m1_w, m2_w)
anova(m2_w, m3_w)

summary(m3_w)  # full model (within-individual slopes by call type)

## 3) Collinearity check 
m_fixed_w <- lm(AOL_RMS ~ NL_c * CALL_TYPE, data = pwdata)
performance::check_collinearity(m_fixed_w)

## 4) Quadratic check
m3_w_quad <- update(m3_w, . ~ . + I(NL_c^2) + I(NL_c^2):CALL_TYPE)
anova(m3_w, m3_w_quad)

## 5) Diagnostics
set.seed(1)
sim_m3_w <- simulateResiduals(m3_w, n = 500, re.form = NULL, plot = TRUE)
par(mfrow = c(1,3))
hist(residuals(m3_w), main = "Residuals histogram", xlab = "Residuals")
qqnorm(residuals(m3_w)); qqline(residuals(m3_w))
plot(residuals(m3_w) ~ fitted(m3_w), pch = 19, col = rgb(0,0,0,0.4),
     xlab = "Fitted", ylab = "Residuals", main = "Residuals vs fitted")
par(mfrow = c(1,1))

## 6) Slopes per call type
ct_within <- emtrends(m3_w_quad, ~ CALL_TYPE, var = "NL_c")
summary(ct_within, infer = TRUE)   

## 7) Overall compensation across call types from the same model
overall_eq <- emtrends(m1_w, ~ 1, var = "NL_c", weights = "equal")
overall_wt <- emtrends(m3_w, ~ 1, var = "NL_c", weights = "proportional")
summary(overall_eq, infer = TRUE)
summary(overall_wt, infer = TRUE)

## 8) Descriptive stats
mf <- model.frame(m3_w)

ct_summary <- mf %>%
  group_by(CALL_TYPE) %>%
  summarise(
    n_calls  = n(),
    mean_AOL = mean(AOL_RMS, na.rm = TRUE),
    sd_AOL   = sd(AOL_RMS,   na.rm = TRUE)
  ) 
ct_summary


