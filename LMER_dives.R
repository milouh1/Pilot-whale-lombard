## Mixed-effects analysis: Dive context
## ==========================================================================

## 0) Load packages
library(lme4)
library(lmerTest)
library(DHARMa)
library(dplyr)
library(emmeans)
library(performance)

emm_options(lmer.df = "satterthwaite")

# 1) Load data 
pw <- read.csv("~/Desktop/stats_PW/Hegeman_Tarifa_PW_Parameters_2025_05_07.csv",
               stringsAsFactors = FALSE, fileEncoding = "UTF-8") %>%
  filter(SNR >= 6) %>%
  group_by(TAG) %>% filter(n() >= 20) %>% ungroup() %>%
  filter(CALL_TYPE == "sk") %>%                     
  mutate(
    TAG = factor(TAG),
    
    # Original contexts
    DIVE_CONTEXT = factor(
      DIVE_CONTEXT,
      levels = c("Surface","Deep_Dive_Descent","Deep_Dive_Ascent","Shallow_Dive")
    ),
    
    # Lump Surface + Shallow_Dive category 
    DIVE_CONTEXT_L = case_when(
      as.character(DIVE_CONTEXT) %in% c("Surface", "Shallow_Dive") ~ "Surface_or_Shallow",
      as.character(DIVE_CONTEXT) == "Deep_Dive_Descent"            ~ "Deep_Dive_Descent",
      as.character(DIVE_CONTEXT) == "Deep_Dive_Ascent"             ~ "Deep_Dive_Ascent",
      TRUE ~ NA_character_
    ),
    DIVE_CONTEXT_L = factor(
      DIVE_CONTEXT_L,
      levels = c("Surface_or_Shallow","Deep_Dive_Descent","Deep_Dive_Ascent")
    ),
    
    # Depth transform 
    AdjustedDepth = log10(ifelse(is.na(DEPTH), NA, pmax(DEPTH, 80)))
  ) %>%
  filter(!is.na(DIVE_CONTEXT_L))

# Check for data 
print(table(pw$DIVE_CONTEXT, useNA = "ifany"))
print(table(pw$DIVE_CONTEXT_L, useNA = "ifany"))

# Within-TAG centering of NL 
pw <- pw %>%
  group_by(TAG) %>%
  mutate(
    NL_bar = mean(NL, na.rm = TRUE), 
    NL_c   = NL - NL_bar
  ) %>%
  ungroup()

# 2) Fitting the models 
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 60000))

m0_w <- lmer(AOL_RMS ~ 1 + (1 + NL_c | TAG),
             data = pw, REML = FALSE, control = ctrl)

m1_w <- lmer(AOL_RMS ~ NL_c + (1 + NL_c | TAG),
             data = pw, REML = FALSE, control = ctrl)

m2_w <- lmer(AOL_RMS ~ NL_c + DIVE_CONTEXT_L + (1 + NL_c | TAG),
             data = pw, REML = FALSE, control = ctrl)

m3a_w <- lmer(AOL_RMS ~ NL_c * DIVE_CONTEXT_L + (1 + NL_c | TAG),
              data = pw, REML = FALSE, control = ctrl)

m3b_w <- lmer(AOL_RMS ~ NL_c + DIVE_CONTEXT_L + AdjustedDepth + (1 + NL_c | TAG),
              data = pw, REML = FALSE, control = ctrl)

## Full model
m4_w <- lmer(AOL_RMS ~ NL_c * DIVE_CONTEXT_L + AdjustedDepth + (1 + NL_c | TAG),
             data = pw, REML = FALSE, control = ctrl)

# Likelihood-ratio tests
anova(m0_w, m1_w)
anova(m1_w, m2_w)
anova(m2_w, m3a_w)
anova(m2_w, m3b_w)
anova(m3b_w, m4_w)
anova(m3a_w, m4_w)

# 3) Collinearity check 
cat("\n--- Collinearity (fixed effects) ---\n")
print(performance::check_collinearity(m4_w))

# 4) Diagnostics 
set.seed(1)
simulateResiduals(m4_w, n = 500, re.form = NULL, plot = TRUE)

par(mfrow = c(1,3))
hist(residuals(m4_w), main = "Residuals histogram", xlab = "Residuals")
qqnorm(residuals(m4_w)); qqline(residuals(m4_w))
plot(residuals(m4_w) ~ fitted(m4_w), pch = 19, col = rgb(0,0,0,0.35),
     xlab = "Fitted", ylab = "Residuals", main = "Residuals vs fitted")
par(mfrow = c(1,1))

# 6) Ns from the data used by m4_w
mf_final <- model.frame(m4_w)

print(
  mf_final %>%
    group_by(DIVE_CONTEXT_L) %>%
    summarise(N_calls = n(),
              N_TAG   = n_distinct(TAG),
              .groups = "drop") %>%
    arrange(desc(N_calls))
)

# 7) Summary of AOL_RMS per dive context
aol_by_ctx <- mf_final %>%
  group_by(DIVE_CONTEXT_L) %>%
  summarise(
    n_calls  = n(),
    n_TAG    = n_distinct(TAG),
    mean_AOL = mean(AOL_RMS, na.rm = TRUE),
    sd_AOL   = sd(AOL_RMS,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_AOL))

print(aol_by_ctx)
