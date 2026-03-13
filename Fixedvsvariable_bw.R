## Comparison of fixed and variable bandwidth: overall compensation to noise 
## ==========================================================================

## 0) Load packages
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 60000))

## 1) Fixed bandwidth
pw_fixedbw <- read.csv("~/Desktop/stats_pw/Hegeman_Tarifa_PW_Parameters_fixedBW_2026_01_10.csv",
                       stringsAsFactors = FALSE)

pw_fixedbw <- pw_fixedbw %>%
  filter(SNR >= 6) %>%                 
  group_by(TAG) %>% filter(n() >= 20) %>% ungroup() %>%
  mutate(
    TAG = factor(TAG),
    CALL_TYPE = factor(CALL_TYPE)
  )

pw_fixedbw <- pw_fixedbw %>%
  group_by(TAG) %>%
  mutate(
    NL_bar = mean(NL, na.rm = TRUE),
    NL_c   = NL - NL_bar
  ) %>%
  ungroup()

m0_fbw_w <- lmer(AOL_RMS ~ 1 + (1 + NL_c | TAG), data = pw_fixedbw, REML = FALSE, control = ctrl)
m1_fbw_w <- lmer(AOL_RMS ~ NL_c + (1 + NL_c | TAG), data = pw_fixedbw, REML = FALSE, control = ctrl)

fixef(m1_fbw_w)["NL_c"]

## 5) Variable bandwidth 

pw_variablebw <- read.csv("~/Desktop/stats_pw/Hegeman_Tarifa_PW_Parameters_2025_05_07.csv",
                          stringsAsFactors = FALSE)

pw_variablebw <- pw_variablebw %>%
  filter(SNR >= 6) %>%
  group_by(TAG) %>% filter(n() >= 20) %>% ungroup() %>%
  mutate(TAG = factor(TAG)) %>%
  group_by(TAG) %>%
  mutate(NL_bar = mean(NL, na.rm = TRUE),
         NL_c   = NL - NL_bar) %>%
  ungroup()

m0_vbw_w <- lmer(AOL_RMS ~ 1 + (1 + NL_c | TAG), data = new, REML = FALSE, control = ctrl)
m1_vbw_w <- lmer(AOL_RMS ~ NL_c + (1 + NL_c | TAG), data = new, REML = FALSE, control = ctrl)

fixef(m1_vbw_w)["NL_c"]

