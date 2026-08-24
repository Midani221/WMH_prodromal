library(dplyr)
library(lubridate)
library(fdapace)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(cluster)
library(lme4)
library(ggeffects)

#This is true
prodromal_data <- "PATH TO FINAL SHEET" # change to sheet


prodromal_data <- read.csv(prodromal_data)

## Descriptive details

prodromal_data <- prodromal_data %>%
  filter(Dx == 4)


prodromal_data$Visit_date_clinical <- mdy(prodromal_data$Visit_date_clinical)

prodromal_data <- prodromal_data %>%
  mutate(
    phenoconversion_date = mdy(phenoconversion_date),
    
    age_at_phenoconversion = if_else(
      !is.na(phenoconversion_date),
      Age + as.numeric(phenoconversion_date - Visit_date_clinical) / 365.25,
      NA_real_
    )
  )

dat <- prodromal_data

subject_summary <- dat %>%
  group_by(SubjID) %>%
  arrange(Age, Visit_date_clinical, .by_group = TRUE) %>%
  summarise(
    n_visits = n_distinct(Visit_No, na.rm = TRUE),
    
    baseline_age = first(Age),
    last_followup_age = last(Age),
    
    baseline_date = if ("Visit_date_clinical" %in% names(dat)) first(Visit_date_clinical) else as.Date(NA),
    last_followup_date = if ("Visit_date_clinical" %in% names(dat)) last(Visit_date_clinical) else as.Date(NA),
    
    followup_years = case_when(
      !is.na(baseline_date) & !is.na(last_followup_date) ~ 
        as.numeric(last_followup_date - baseline_date) / 365.25,
      TRUE ~ last_followup_age - baseline_age
    ),
    
    baseline_UPDRS3 = first(MDS_UPDRS3_OFF_SCORE),
    baseline_scopa = first(scopa),
    baseline_MoCA = first(MoCA),
    baseline_Cognition_EF = first(Cognition_EF),
    baseline_Cognition_M = first(Cognition_M),
    
    phenoconverted = first(phenoconversion),
    age_at_phenoconversion  = first(age_at_phenoconversion),
    .groups = "drop"
  )


median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  med <- median(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  
  paste0(
    round(med, digits),
    " [",
    round(q1, digits),
    ", ",
    round(q3, digits),
    "]"
  )
}


descriptive_table <- tibble(
  Characteristic = c(
    "Number of participants",
    "Number of visits",
    "Median visits per participant",
    "Median follow-up time, years",
    "Age at baseline, years",
    "Age at last follow-up, years",
    "Baseline MDS-UPDRS III",
    "Baseline SCOPA-AUT",
    "Baseline MoCA",
    "Baseline executive function composite",
    "Baseline memory composite",
    "Number converting to PD",
    "Percentage converting to PD",
    "Age at phenoconversion, years"
  ),
  
  Value = c(
    n_distinct(subject_summary$SubjID),
    
    sum(subject_summary$n_visits, na.rm = TRUE),
    
    median_iqr(subject_summary$n_visits, digits = 0),
    
    median_iqr(subject_summary$followup_years),
    
    median_iqr(subject_summary$baseline_age),
    
    median_iqr(subject_summary$last_followup_age),
    
    median_iqr(subject_summary$baseline_UPDRS3),
    
    median_iqr(subject_summary$baseline_scopa),
    
    median_iqr(subject_summary$baseline_MoCA),
    
    median_iqr(subject_summary$baseline_Cognition_EF),
    
    median_iqr(subject_summary$baseline_Cognition_M),
    
    sum(subject_summary$phenoconverted, na.rm = TRUE),
    
    paste0(
      round(
        100 * mean(subject_summary$phenoconverted, na.rm = TRUE),
        1
      ),
      "%"
    ),
    
    median_iqr(subject_summary$age_at_phenoconversion)
  )
)

descriptive_table



## Quick Spaghetti Plots for MDRS3, SCOPA, MoCA, EF

prodromal_updrs <- prodromal_data %>%
  filter(Dx == 4) %>%
  filter(!is.na(SubjID),
         !is.na(Age)) %>%
  mutate(
    phenoconversion_group = ifelse(phenoconversion == 1,
                                   "Phenoconverter",
                                   "Non-converter")
  ) %>%
  arrange(SubjID, Age)



prodromal_updrs$WMH_total <- as.numeric(prodromal_updrs$WMH_total_lin*1000*100/(prodromal_updrs$scaling_factor*prodromal_updrs$Native_v))
prodromal_updrs$WMH_total <- log(1+ prodromal_updrs$WMH_total)

prodromal_updrs$dWMH <- as.numeric(prodromal_updrs$WMH_deep_lin*1000*100/(prodromal_updrs$scaling_factor*prodromal_updrs$Native_v))
prodromal_updrs$dWMH <- log(1 + as.numeric(prodromal_updrs$dWMH))

prodromal_updrs$pWMH <- as.numeric(prodromal_updrs$WMH_peri_lin*1000*100/(prodromal_updrs$scaling_factor*prodromal_updrs$Native_v))
prodromal_updrs$pWMH<- log(1+ as.numeric(prodromal_updrs$pWMH))


#prodromal_updrs <- prodromal_updrs %>%
#filter(QC_final == TRUE)

vars <- c(
  "MDS_UPDRS3_OFF_SCORE",
  "scopa",
  "MoCA",
  "Cognition_EF",
  "Cognition_M",
  "WMH_total"
)

var_labels <- c(
  MDS_UPDRS3_OFF_SCORE = "MDS-UPDRS III OFF score",
  scopa = "SCOPA-AUT score",
  MoCA = "MoCA score",
  Cognition_EF = "Executive function composite",
  Cognition_M = "Memory composite",
  WMH_total = "Total WMHs"
)

make_spaghetti_plot <- function(data, outcome_var) {
  
  y_label <- var_labels[[outcome_var]]
  if (is.null(y_label)) y_label <- outcome_var
  
  data_plot <- data %>%
    filter(
      !is.na(years_from_baseline),
      !is.na(.data[[outcome_var]]),
      !is.na(SubjID)    )
  
  ggplot(
    data_plot,
    aes(
      x = years_from_baseline,
      y = .data[[outcome_var]],
      group = SubjID
      ,color = phenoconversion_group
    )
  ) +
    geom_line(alpha = 0.25, linewidth = 0.4) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(
      aes(group = phenoconversion_group),
      method = "loess",
     se = TRUE,
      linewidth = 1.2
    ) +
    labs(
      title = paste0(y_label, " trajectories by years from baseline in prodromal participants"),
      x = "Years from baseline",
      y = y_label
      ,color = "Phenoconversion status"
    ) +
    theme_classic()
}



#prodromal_pheno <- prodromal_updrs[prodromal_updrs$phenoconversion_group == "Phenoconverter",]

#prodromal_pheno$age_p <- prodromal_pheno$age_at_phenoconversion - prodromal_pheno$Age



spaghetti_plots <- map(vars, ~ make_spaghetti_plot(prodromal_model_data, .x))

names(spaghetti_plots) <- vars

spaghetti_plots$MDS_UPDRS3_OFF_SCORE
spaghetti_plots$scopa
spaghetti_plots$MoCA
spaghetti_plots$Cognition_EF
spaghetti_plots$Cognition_M
spaghetti_plots$WMH_total


#### Fix data a bit - start here
prodromal_long <- prodromal_data %>%  
  mutate(
    Age = as.numeric(Age),
    scopa = as.numeric(as.character(scopa))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(Age)
  ) %>%
  arrange(SubjID, Age)


## Baseline data
# 2. Create one-row-per-subject baseline covariate table
baseline_covariates <- prodromal_long %>%
  group_by(SubjID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    subgroup_clean = str_squish(coalesce(Subgroup, "")),

    RBD      = as.factor(str_detect(subgroup_clean, "\\bRBD\\b")),
    LRRK2    = as.factor(str_detect(subgroup_clean, "\\bLRRK2\\b")),
    GBA      = as.factor(str_detect(subgroup_clean, "\\bGBA\\b")),
    SNCA     = as.factor(str_detect(subgroup_clean, "\\bSNCA\\b")),
    PRKN     = as.factor(str_detect(subgroup_clean, "\\bPRKN\\b")),
    Hyposmia = as.factor(str_detect(subgroup_clean, "\\bHyposmia\\b")),

    baseline_age = Age,
    baseline_scopa = scopa,

    phenoconversion_group = ifelse(
      phenoconversion == 1,
      "Phenoconverter",
      "Non-converter"
    )
  ) %>%
  select(
    SubjID,
    baseline_age,
    baseline_scopa,
    subgroup_clean,
    RBD,
    LRRK2,
    GBA,
    SNCA,
    PRKN,
    Hyposmia,
    phenoconversion_group,

  )


prodromal_model_data <- prodromal_long %>%
  left_join(baseline_covariates, by = "SubjID")


prodromal_model_data <- prodromal_model_data %>%
  mutate(
    Age = as.numeric(Age)
  ) %>%
  arrange(SubjID, Age) %>%
  group_by(SubjID) %>%
  mutate(
    baseline_age = first(Age),
    years_from_baseline = Age - baseline_age
  ) %>%
  ungroup() %>%
  mutate(
    baseline_age_c = baseline_age - mean(baseline_age, na.rm = TRUE)
  )

prodromal_model_data$Native_v <- prodromal_model_data$Native_v/1000

prodromal_model_data$WMH_total <- as.numeric(prodromal_model_data$WMH_native_roi)
prodromal_model_data$WMH_total <- log(1+ prodromal_model_data$WMH_total)


prodromal_model_data$dWMH <- as.numeric(prodromal_model_data$WMH_deep_roi)
prodromal_model_data$dWMH <- log(1 + as.numeric(prodromal_model_data$dWMH))

prodromal_model_data$pWMH <- as.numeric(prodromal_model_data$WMH_peri_roi)
prodromal_model_data$pWMH<- log(1+ as.numeric(prodromal_model_data$pWMH))

hist(prodromal_model_data$Native_v)

prodromal_model_data <- prodromal_model_data %>%
   filter(Age >= 50, Age <= 85) %>%
   group_by(SubjID) %>%
   ungroup()


# Variables: MDS_UPDRS3_OFF_SCORE, MoCA, scopa, hyposmia, orthostasis, Cognition_EF, Cognition_M, 
# Cognition_L, Cognition_VF, Cognition_AWM

updrs_model <- lmer(
  scopa ~ Age + Sx + Education_Years + phenoconversion_group*pWMH + Native_v + Equivalent_Vascular_RF +
    (1 | SubjID),
  data = prodromal_model_data,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

#model <- lmer(
#  WMH_total ~ Age + (1 | SubjID),
#  data = prodromal_model_data,
#  REML = FALSE
#)

#summary(model)

summary(updrs_model)

coef(summary(updrs_model))

confint(updrs_model, method = "Wald")


###
#WMH model


wmh_model <- lmer(
  WMH_total ~
    Age*phenoconversion_group +
     Native_v +
    Sx +
    (1 | SubjID),
  data = prodromal_model_data,
  REML = FALSE
)


summary(wmh_model)

confint(wmh_model)

## AGE at visit analysis:



# -----------------------------
# 1. Clean data
# -----------------------------

prodromal_model_data$hyposmia <- as.factor(prodromal_model_data$hyposmia)
prodromal_model_data$orthostasis <- as.factor(prodromal_model_data$orthostasis)



baseline_chars <- prodromal_model_data %>%
  select(SubjID, years_from_baseline, Age, Sx, scopa, phenoconversion_group, Equivalent_Vascular_RF, hyposmia, orthostasis, Education_Years, WMH_total,dWMH,pWMH,Native_v) %>%
  mutate(
    Age = as.numeric(Age),
    baseline_scopa = as.numeric(as.character(scopa))
  ) %>%
  filter(!is.na(SubjID), !is.na(Age)) %>%
  arrange(SubjID, Age) %>%
  group_by(SubjID) %>%
  slice(1) %>%
  ungroup() %>%
  rename(
    baseline_age = Age
  ) 

# FPC plot function


summary(as.factor(round(prodromal_model_data$years_from_baseline)))


## Prepare FPC data
dat <- prodromal_model_data %>%
  select(SubjID, years_from_baseline, MDS_UPDRS3_OFF_SCORE) %>%
  mutate(
    years_from_baseline = as.numeric(years_from_baseline),
    UPDRS3 = as.numeric(as.character(MDS_UPDRS3_OFF_SCORE))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(years_from_baseline),
    !is.na(UPDRS3)
  ) %>%
  arrange(SubjID, years_from_baseline)

# -----------------------------
# 2. Resolve duplicate visits at same age
# -----------------------------
dat <- dat %>%
  group_by(SubjID, years_from_baseline) %>%
  summarise(
    UPDRS3 = mean(UPDRS3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SubjID, years_from_baseline)

# -----------------------------
# 3. Choose analysis dataset
# -----------------------------

dat2 <- dat %>%
 filter(years_from_baseline < 10.5)

# -----------------------------
# 4. Create FPCA lists safely
# -----------------------------
fpca_lists <- dat2 %>%
  arrange(SubjID, years_from_baseline) %>%
  group_by(SubjID) %>%
  summarise(
    Lt = list(years_from_baseline),
    Ly = list(UPDRS3),
    n_visits = n(),
    .groups = "drop"
  )

ids <- fpca_lists$SubjID
Lt  <- fpca_lists$Lt
Ly  <- fpca_lists$Ly

# -----------------------------
# 5. Check inputs
# -----------------------------
stopifnot(length(ids) == length(Lt))
stopifnot(length(ids) == length(Ly))
stopifnot(all(lengths(Lt) == lengths(Ly)))

summary(lengths(Lt))
range(unlist(Lt), na.rm = TRUE)
range(unlist(Ly), na.rm = TRUE)

# -----------------------------
# 6. Run FPCA
# -----------------------------
optns <- list(
  dataType = "Sparse",
  kernel = "gauss",
  error = TRUE,
  methodSelectK = "FVE",
  FVEthreshold = 0.95,
  maxK = 5,
  verbose = TRUE
)

fpca_fit <- FPCA(Ly = Ly, Lt = Lt, optns = optns)

CreateDesignPlot(Lt)
plot(fpca_fit)

# -----------------------------
# 7. Extract first two scores
# -----------------------------

## Age at visit analysis looks great.

scores <- as.data.frame(fpca_fit$xiEst)
colnames(scores) <- paste0("FPC", seq_len(ncol(scores)))
scores$SubjID <- ids

head(scores)

### Plotting graph and getting data for patients

grid <- fpca_fit$workGrid
mu   <- fpca_fit$mu
phi  <- fpca_fit$phi[, 1:3, drop = FALSE]


fpca_merged <- scores %>%
  left_join(baseline_chars, by = "SubjID")

model_fpc1 <- lm(FPC1 ~ baseline_age + Sx + baseline_scopa + hyposmia + orthostasis + WMH_total,
                 data = fpca_merged)


summary(model_fpc1)

model_fpc2 <- lm(FPC2 ~  baseline_age + Sx + baseline_scopa + hyposmia + orthostasis + WMH_total,
                 data = fpca_merged)

summary(model_fpc2)





fpca_merged$phenoconversion_group <- as.factor(fpca_merged$phenoconversion_group)


fpca_merged <- fpca_merged %>%
  mutate(
    phenoconversion_binary = case_when(
      phenoconversion_group %in% c("Converter", "Phenoconverter", "1", 1) ~ 1,
      phenoconversion_group %in% c("Non-converter", "Nonconverter", "0", 0) ~ 0,
      TRUE ~ NA_real_
    )
  )

model_fpc1_x <- lm(FPC1 ~ baseline_age + Sx + baseline_scopa + RBD + LRRK2 + GBA + Hyposmia + phenoconversion_group + log_deep_wmh_ml +log_periventricular_wmh_ml,
                  data = fpca_merged3)

summary(model_fpc1_x)
lambda <- fpca_fit$lambda

fve_table <- data.frame(
  FPC = paste0("FPC", seq_along(lambda)),
  Eigenvalue = lambda,
  FVE = lambda / sum(lambda),
  Cumulative_FVE = cumsum(lambda / sum(lambda))
)

fve_table

# Save the fpca_fit for subsequent analysis
saveRDS(fpca_fit, "path
")

# Plot FPCs to investigate how they affect the data
plot_fpc_effect <- function(k) {
  
  lambda_k <- fpca_fit$lambda[k]
  phi_k <- fpca_fit$phi[, k]
  
  plot_df <- data.frame(
    age = fpca_fit$workGrid,
    mean = fpca_fit$mu,
    plus = fpca_fit$mu + 2 * sqrt(lambda_k) * phi_k,
    minus = fpca_fit$mu - 2 * sqrt(lambda_k) * phi_k
  )
  
  ggplot(plot_df, aes(x = age)) +
    geom_line(aes(y = mean), linewidth = 1.2) +
    geom_line(aes(y = plus), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = minus), linetype = "dashed", linewidth = 1) +
    labs(
      title = paste0("Interpretation of FPC", k),
      x = "Age at visit",
      y = "MDS-UPDRS III OFF score"
    ) +
    theme_classic()
}

plot_fpc_effect(1)
plot_fpc_effect(2)
plot_fpc_effect(3)

#
ggplot(fpca_merged, aes(x = phenoconversion_group, y = WMH_total)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  labs(
    title = "Log WMH by phenoconversion status",
    x = "Phenoconversion status",
    y = "Log WMH percentage"
  ) +
  theme_classic()

ggplot(fpca_merged, aes(x = phenoconversion_group, y = FPC2)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  labs(
    title = "FPC2 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "FPC2 score"
  ) +
  theme_classic()

ggplot(fpca_merged, aes(x = phenoconversion_group, y = FPC3)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  labs(
    title = "FPC3 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "FPC3 score"
  ) +
  theme_classic()


ggplot(
  fpca_merged3,
  aes(x = FPC1, y = log_periventricular_wmh_ml, color = phenoconversion_group)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "UPDRS3 FPCA score space",
    x = "FPC1",
    y = "FPC2",
    color = "Phenoconversion status"
  ) +
  theme_classic()
###



fpca_merged <- fpca_merged %>%
  mutate(
    FPC1_z = as.numeric(scale(FPC1)),
    FPC2_z = as.numeric(scale(FPC2))
  )

fpca_merged$phenoconversion_group <- ifelse(fpca_merged$phenoconversion_group == "Phenoconverter", 1, 0)

fpca_merged$phenoconversion_group <- ifelse(fpca_merged$phenoconversion_group == 1, "Phenoconverter", "Non-converter")

fpca_merged$phenoconversion_group <- relevel(
  factor(fpca_merged$phenoconversion_group),
  ref = "Non-converter"
)


fpca_merged <- as.data.frame(fpca_merged)
model_pheno_z <- glm(
  phenoconversion_group ~  FPC1_z + FPC2_z + baseline_age + Sx + orthostasis + hyposmia + baseline_scopa + WMH_total ,
  data = fpca_merged,
  family = binomial
)

summary(model_pheno_z)



### SCOPA
dat <- prodromal_data %>%
  select(SubjID, years_from_baseline, scopa) %>%
  mutate(
    age_at_visit = as.numeric(Age),
    SCOPA = as.numeric(as.character(scopa))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(age_at_visit),
    !is.na(SCOPA)
  ) %>%
  arrange(SubjID, age_at_visit)

# -----------------------------
# 2. Resolve duplicate visits at same age
# -----------------------------
dat <- dat %>%
  group_by(SubjID, age_at_visit) %>%
  summarise(
    SCOPA = mean(SCOPA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SubjID, age_at_visit)

# -----------------------------
# 3. Choose analysis dataset
# -----------------------------

# Option B: restricted age range
 dat2 <- dat %>%
   filter(age_at_visit >= 50, age_at_visit <= 85) %>%
   group_by(SubjID) %>%
   ungroup()

# -----------------------------
# 4. Create FPCA lists safely
# -----------------------------
fpca_lists <- dat2 %>%
  arrange(SubjID, age_at_visit) %>%
  group_by(SubjID) %>%
  summarise(
    Lt = list(age_at_visit),
    Ly = list(SCOPA),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  filter(n_visits >= 2)

ids <- fpca_lists$SubjID
Lt  <- fpca_lists$Lt
Ly  <- fpca_lists$Ly

# -----------------------------
# 5. Check inputs
# -----------------------------
stopifnot(length(ids) == length(Lt))
stopifnot(length(ids) == length(Ly))
stopifnot(all(lengths(Lt) == lengths(Ly)))

summary(lengths(Lt))
range(unlist(Lt), na.rm = TRUE)
range(unlist(Ly), na.rm = TRUE)

# -----------------------------
# 6. Run FPCA
# -----------------------------
optns <- list(
  dataType = "Sparse",
  kernel = "gauss",
  error = TRUE,
  methodSelectK = "FVE",
  FVEthreshold = 0.95,
  maxK = 5,
  verbose = TRUE
)

fpca_fit_scopa <- FPCA(Ly = Ly, Lt = Lt, optns = optns)
fpca_fit_scopa <- readRDS("path")


plot(fpca_fit_scopa)

scores_scopa <- as.data.frame(fpca_fit_scopa$xiEst)
colnames(scores_scopa) <- paste0("FPC_scopa", seq_len(ncol(scores_scopa)))
scores_scopa$SubjID <- ids


### Plotting graph and getting data for patients

grid <- fpca_fit_scopa$workGrid
mu   <- fpca_fit_scopa$mu
phi  <- fpca_fit_scopa$phi[, 1:3, drop = FALSE]


fpca_merged2 <- fpca_merged %>%
  left_join(scores_scopa, by = "SubjID")

model_fpc1 <- lm(FPC_scopa1_z ~ baseline_age + Sx  + orthostasis + hyposmia + WMH_total ,
                 data = fpca_merged2)


summary(model_fpc1)

model_fpc2 <- lm(FPC_scopa2 ~  baseline_age + Sx  + hyposmia + orthostasis + WMH_total,
                 data = fpca_merged2)

summary(model_fpc2)


model_fpc3 <- lm(FPC_scopa3 ~  baseline_age + Sx  + RBD + LRRK2 + GBA + Hyposmia + phenoconversion_group,
                 data = fpca_merged2)


summary(model_fpc3)

## Plot scopa
plot_fpc_effect <- function(k,fpca) {
  fpca <- fpca
  lambda_k <- fpca$lambda[k]
  phi_k <- fpca$phi[, k]
  
  plot_df <- data.frame(
    age = fpca$workGrid,
    mean = fpca$mu,
    plus = fpca$mu + 2 * sqrt(lambda_k) * phi_k,
    minus = fpca$mu - 2 * sqrt(lambda_k) * phi_k
  )
  
  ggplot(plot_df, aes(x = age)) +
    geom_line(aes(y = mean), linewidth = 1.2) +
    geom_line(aes(y = plus), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = minus), linetype = "dashed", linewidth = 1) +
    labs(
      title = paste0("Interpretation of FPC", k),
      x = "Age at visit",
      y = "Variable Score"
    ) +
    theme_classic()
}

plot_fpc_effect(1,fpca_fit_scopa)
plot_fpc_effect(2,fpca_fit_scopa)
plot_fpc_effect(3, fpca_fit_scopa)


fpca_merged2 <- fpca_merged2 %>%
  mutate(
    FPC_scopa1_z = as.numeric(scale(FPC_scopa1)),
    FPC_scopa2_z = as.numeric(scale(FPC_scopa2)),
    FPC_scopa3_z = as.numeric(scale(FPC_scopa3))
  )

#
ggplot(fpca_merged2, aes(x = phenoconversion_group, y = FPC_scopa2)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  labs(
    title = "FPC_scopa1 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "FPC_scopa1 score"
  ) +
  theme_classic()

ggplot(fpca_merged2, aes(x = phenoconversion_group, y = FPC2)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  labs(
    title = "FPC2 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "FPC2 score"
  ) +
  theme_classic()




ggplot(
  fpca_merged2,
  aes(x = FPC1, y = FPC_scopa1, color = phenoconversion_group)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "UPDRS3 FPCA score space",
    x = "FPC1",
    y = "FPC_scopa1",
    color = "Phenoconversion status"
  ) +
  theme_classic()


model_fpcx1 <- lm(FPC1 ~  FPC_scopa1 + FPC_scopa2 + baseline_age + Sx + baseline_scopa + RBD + LRRK2 + GBA + Hyposmia + phenoconversion_group,
                 data = fpca_merged2)

summary(model_fpcx1)

model_fpcx2 <- lm(FPC2 ~  FPC_scopa1 + FPC_scopa2 + baseline_age + Sx + baseline_scopa + RBD + LRRK2 + GBA + Hyposmia + phenoconversion_group,
                 data = fpca_merged2)

summary(model_fpcx2)

saveRDS(fpca_fit_scopa, "path")


##

fpca_merged2$phenoconversion_group <- ifelse(fpca_merged2$phenoconversion_group == "Phenoconverter", 1, 0)

fpca_merged2$phenoconversion_group <- ifelse(fpca_merged2$phenoconversion_group == 1, "Phenoconverter", "Non-converter")

model_glm <- glm(phenoconversion_group ~ FPC1 + FPC2 + FPC_scopa1 + FPC_scopa2 + 
baseline_age + Sx + WMH_total, family = "binomial",
data = fpca_merged2)
summary(model_glm)
## Analysis of FPCs EFFECT


##

prodromal_updrs <- prodromal_data %>%
  filter(Dx == 4) %>%
  filter(!is.na(SubjID),
         !is.na(Age),
         !is.na(MDS_UPDRS3_OFF_SCORE)) %>%
  mutate(
    phenoconversion_group = ifelse(phenoconversion == 1,
                                   "Phenoconverter",
                                   "Non-converter")
  ) %>%
  arrange(SubjID, Age)

ggplot(prodromal_updrs,
       aes(x = Age,
           y = MDS_UPDRS3_OFF_SCORE,
           group = SubjID,
           color = phenoconversion_group)) +
  geom_line(alpha = 0.25, linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(
    aes(group = phenoconversion_group),
    method = "loess",
    se = TRUE,
    linewidth = 1.2
  ) +
  labs(
    title = "UPDRS III trajectories by age in prodromal participants",
    x = "Age at visit",
    y = "UPDRS III score",
    color = "Phenoconversion status"
  ) +
  theme_classic()

## Plot for one subject

# -----------------------------
# Identify subject in FPCA object
# -----------------------------

target_id <- 75503

ids <- as.character(ids)

subj_index <- match(target_id, ids)

if (is.na(subj_index)) {
  stop(
    paste0(
      "Subject ", target_id, " was not found in ids. ",
      "Check unique(ids) to confirm the exact subject ID format."
    )
  )
}

# -----------------------------
# Choose how many FPCs to use
# -----------------------------

# Use all selected FPCA components
K <- ncol(fpca_fit$xiEst)

# Or, if you only want first 3 components, use:
# K <- min(3, ncol(fpca_fit$xiEst))

# -----------------------------
# Extract FPCA objects
# -----------------------------

grid <- fpca_fit$workGrid
mu   <- fpca_fit$mu
phi  <- fpca_fit$phi[, 1:K, drop = FALSE]
xi   <- fpca_fit$xiEst[subj_index, 1:K]

# -----------------------------
# Reconstruct subject-specific trajectory
# -----------------------------

subject_fit <- as.numeric(mu + phi %*% xi)

fpca_plot_df <- tibble(
  age_at_visit = grid,
  population_mean = as.numeric(mu),
  subject_fpca_fit = subject_fit
)

# -----------------------------
# Extract observed data for this subject
# -----------------------------

obs_df <- dat2 %>%
  mutate(SubjID = as.character(SubjID)) %>%
  filter(SubjID == target_id) %>%
  arrange(age_at_visit)

if (nrow(obs_df) == 0) {
  stop(
    paste0(
      "No observed data found for subject ", target_id, 
      " in dat2. Check whether the subject was removed by the age filter."
    )
  )
}

# -----------------------------
# Plot
# -----------------------------

ggplot() +
  geom_line(
    data = fpca_plot_df,
    aes(x = age_at_visit, y = population_mean),
    linetype = "dashed",
    linewidth = 1,
    alpha = 0.7
  ) +
  geom_line(
    data = fpca_plot_df,
    aes(x = age_at_visit, y = subject_fpca_fit),
    linewidth = 1.2
  ) +
  geom_point(
    data = obs_df,
    aes(x = age_at_visit, y = UPDRS3),
    size = 3
  ) +
  geom_line(
    data = obs_df,
    aes(x = age_at_visit, y = UPDRS3),
    linewidth = 0.7,
    alpha = 0.5
  ) +
  labs(
    title = paste("MDS-UPDRS III OFF trajectory for subject", target_id),
    subtitle = paste0(

      "Observed scores and FPCA-reconstructed trajectory using ",
      K, " FPC(s)"
    ),
    x = "Age at visit",
    y = "MDS-UPDRS III OFF score",
    caption = "Solid line = subject-specific FPCA trajectory; dashed line = population mean; points = observed visits"
  ) +
  theme_classic(base_size = 14)

#########3 MoCA
# ============================================================
# 1. Prepare longitudinal MoCA data
# ============================================================

dat <- prodromal_data %>%
  select(SubjID, Age, MoCA) %>%
  mutate(
    age_at_visit = as.numeric(Age),
    MOCA = as.numeric(as.character(MoCA))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(age_at_visit),
    !is.na(MOCA)
  ) %>%
  arrange(SubjID, age_at_visit)

# Resolve duplicate observations at the same age
dat <- dat %>%
  group_by(SubjID, age_at_visit) %>%
  summarise(
    MOCA = mean(MOCA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SubjID, age_at_visit)

# Derive baseline MoCA from the earliest available observation
baseline_moca_df <- dat %>%
  group_by(SubjID) %>%
  slice_min(
    order_by = age_at_visit,
    n = 1,
    with_ties = FALSE
  ) %>%
  transmute(
    SubjID,
    baseline_moca = MOCA
  ) %>%
  ungroup()

# ============================================================
# 2. Restrict the analysis age range
# ============================================================

dat2 <- dat %>%
  filter(
    age_at_visit >= 50,
    age_at_visit <= 85
  )

# ============================================================
# 3. Create FPCA input lists
# ============================================================

fpca_lists <- dat2 %>%
  arrange(SubjID, age_at_visit) %>%
  group_by(SubjID) %>%
  summarise(
    Lt = list(age_at_visit),
    Ly = list(MOCA),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  filter(n_visits >= 2)

ids <- fpca_lists$SubjID
Lt <- fpca_lists$Lt
Ly <- fpca_lists$Ly

stopifnot(length(ids) == length(Lt))
stopifnot(length(ids) == length(Ly))
stopifnot(all(lengths(Lt) == lengths(Ly)))

summary(lengths(Lt))
range(unlist(Lt), na.rm = TRUE)
range(unlist(Ly), na.rm = TRUE)

# ============================================================
# 4. Run or load MoCA FPCA
# ============================================================

optns <- list(
  dataType = "Sparse",
  kernel = "gauss",
  error = TRUE,
  methodSelectK = "FVE",
  FVEthreshold = 0.95,
  maxK = 5,
  verbose = TRUE
)

fpca_file <- paste0(
  "path",
  "fpca_fit_moca.rds"
)

if (file.exists(fpca_file)) {
  fpca_fit_moca <- readRDS(fpca_file)
} else {
  fpca_fit_moca <- FPCA(
    Ly = Ly,
    Lt = Lt,
    optns = optns
  )

  saveRDS(fpca_fit_moca, fpca_file)
}

plot(fpca_fit_moca)

# This code below assumes at least three FPCs were retained
stopifnot(ncol(fpca_fit_moca$xiEst) >= 3)

# ============================================================
# 5. Extract and standardize MoCA FPC scores
# ============================================================

scores_moca <- as.data.frame(fpca_fit_moca$xiEst)

colnames(scores_moca) <- paste0(
  "FPC_moca",
  seq_len(ncol(scores_moca))
)

scores_moca$SubjID <- ids

fpca_merged3 <- fpca_merged2 %>%
  select(-any_of("baseline_moca")) %>%
  left_join(scores_moca, by = "SubjID") %>%
  left_join(baseline_moca_df, by = "SubjID") %>%
  mutate(
    FPC_moca1_z = as.numeric(scale(FPC_moca1)),
    FPC_moca2_z = as.numeric(scale(FPC_moca2)),
    FPC_moca3_z = as.numeric(scale(FPC_moca3))
  )

# ============================================================
# 6. Predictors of MoCA FPC scores
# ============================================================

model_moca_fpc1 <- lm(
  FPC_moca1_z ~
    baseline_age +
    Sx +
    baseline_moca +
    orthostasis +
    hyposmia + baseline_scopa +
    pWMH + dWMH +  phenoconversion_group,
  data = fpca_merged3
)

summary(model_moca_fpc1)

model_moca_fpc2 <- lm(
  FPC_moca2_z ~
    baseline_age +
    Sx +
    baseline_moca +
    orthostasis +
    hyposmia +
    WMH_total,
  data = fpca_merged
)

summary(model_moca_fpc2)

model_moca_fpc3 <- lm(
  FPC_moca1_z ~
    baseline_age +
    Sx +
    baseline_moca +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group +
    WMH_total,
  data = fpca_merged3
)

summary(model_moca_fpc3)

# ============================================================
# 7. Plot the MoCA FPCA components
# ============================================================

plot_fpc_effect <- function(k, fpca, outcome_label = "MoCA score") {

  lambda_k <- fpca$lambda[k]
  phi_k <- fpca$phi[, k]

  plot_df <- data.frame(
    age = fpca$workGrid,
    mean = fpca$mu,
    plus = fpca$mu + 2 * sqrt(lambda_k) * phi_k,
    minus = fpca$mu - 2 * sqrt(lambda_k) * phi_k
  )

  ggplot(plot_df, aes(x = age)) +
    geom_line(
      aes(y = mean),
      linewidth = 1.2
    ) +
    geom_line(
      aes(y = plus),
      linetype = "dotted",
      linewidth = 1
    ) +
    geom_line(
      aes(y = minus),
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste0("Interpretation of MoCA FPC", k),
      subtitle = "Mean trajectory and mean ± 2 SD component variation",
      x = "Age at visit",
      y = outcome_label
    ) +
    theme_classic()
}

plot_fpc_effect(1, fpca_fit_moca)
plot_fpc_effect(2, fpca_fit_moca)
plot_fpc_effect(3, fpca_fit_moca)

# ============================================================
# 8. Compare MoCA FPC scores by phenoconversion
# ============================================================

ggplot(
  fpca_merged3,
  aes(
    x = phenoconversion_group,
    y = FPC_moca1_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "MoCA FPC1 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized MoCA FPC1 score"
  ) +
  theme_classic()

ggplot(
  fpca_merged2,
  aes(
    x = phenoconversion_group,
    y = FPC_moca2_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "MoCA FPC2 scores by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized MoCA FPC2 score"
  ) +
  theme_classic()

# ============================================================
# 9. Relationship between UPDRS and MoCA FPCA scores
# ============================================================

ggplot(
  fpca_merged3,
  aes(
    x = FPC1_z,
    y = FPC_moca1_z,
    color = phenoconversion_group
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  labs(
    title = "Relationship between UPDRS and MoCA FPCA scores",
    x = "UPDRS FPC1 score",
    y = "Standardized MoCA FPC1 score",
    color = "Phenoconversion status"
  ) +
  theme_classic()

# ============================================================
# 10. Cross-modal regression models
# ============================================================

model_fpcx1 <- lm(
  FPC1_z ~
    FPC_moca1_z +
    FPC_moca2_z +
    FPC_scopa1_z +
    FPC_scopa2_z +
    baseline_age +
    Sx +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group +
    WMH_total,
  data = fpca_merged3
)

summary(model_fpcx1)

model_fpcx2 <- lm(
  FPC2 ~
    FPC_moca1_z +
    FPC_moca2_z +
    baseline_age +
    Sx +
    baseline_moca +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group,
  data = fpca_merged2
)

summary(model_fpcx2)

# ============================================================
# 11. Logistic regression for phenoconversion
# ============================================================

fpca_merged2 <- fpca_merged2 %>%
  mutate(
    phenoconversion_group = factor(
      phenoconversion_group,
      levels = c(
        "Non-converter",
        "Phenoconverter"
      )
    )
  )

table(
  fpca_merged2$phenoconversion_group,
  useNA = "ifany"
)

model_glm <- glm(
  phenoconversion_group ~
    FPC1_z +
    FPC2_z +
    FPC_moca1_z +
    FPC_moca2_z +
    FPC_scopa1_z +
    FPC_scopa2_z +
    baseline_age +
    Sx +
    WMH_total,
  family = binomial(link = "logit"),
  data = fpca_merged3
)

summary(model_glm)

# Odds ratios and 95% confidence intervals
glm_results <- data.frame(
  variable = names(coef(model_glm)),
  OR = exp(coef(model_glm)),
  CI_low = exp(
    coef(model_glm) -
      1.96 * summary(model_glm)$coefficients[, "Std. Error"]
  ),
  CI_high = exp(
    coef(model_glm) +
      1.96 * summary(model_glm)$coefficients[, "Std. Error"]
  ),
  p_value = summary(model_glm)$coefficients[, "Pr(>|z|)"]
)

glm_results


###########################

# ============================================================
# 1. Prepare executive-function data
# ============================================================

dat_ef <- prodromal_data %>%
  select(SubjID, Age, Cognition_EF) %>%
  mutate(
    age_at_visit = as.numeric(Age),
    EF = as.numeric(as.character(Cognition_EF))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(age_at_visit),
    !is.na(EF)
  ) %>%
  arrange(SubjID, age_at_visit)

# ============================================================
# 2. Resolve duplicate visits at the same age
# ============================================================

dat_ef <- dat_ef %>%
  group_by(SubjID, age_at_visit) %>%
  summarise(
    EF = mean(EF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SubjID, age_at_visit)

# Earliest available EF value for each participant
baseline_ef <- dat_ef %>%
  group_by(SubjID) %>%
  slice_min(
    order_by = age_at_visit,
    n = 1,
    with_ties = FALSE
  ) %>%
  transmute(
    SubjID,
    baseline_EF = EF
  ) %>%
  ungroup()

# ============================================================
# 3. Restrict age range
# ============================================================

dat2_ef <- dat_ef %>%
  filter(
    age_at_visit >= 50,
    age_at_visit <= 85
  )

# ============================================================
# 4. Create FPCA lists
# ============================================================

fpca_lists_ef <- dat2_ef %>%
  arrange(SubjID, age_at_visit) %>%
  group_by(SubjID) %>%
  summarise(
    Lt = list(age_at_visit),
    Ly = list(EF),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  filter(n_visits >= 2)

ids_ef <- fpca_lists_ef$SubjID
Lt_ef  <- fpca_lists_ef$Lt
Ly_ef  <- fpca_lists_ef$Ly

# ============================================================
# 5. Check inputs
# ============================================================

stopifnot(length(ids_ef) == length(Lt_ef))
stopifnot(length(ids_ef) == length(Ly_ef))
stopifnot(all(lengths(Lt_ef) == lengths(Ly_ef)))

summary(lengths(Lt_ef))
range(unlist(Lt_ef), na.rm = TRUE)
range(unlist(Ly_ef), na.rm = TRUE)

# ============================================================
# 6. Run FPCA
# ============================================================

optns <- list(
  dataType = "Sparse",
  kernel = "gauss",
  error = TRUE,
  methodSelectK = "FVE",
  FVEthreshold = 0.95,
  maxK = 5,
  verbose = TRUE
)

fpca_fit_ef <- FPCA(
  Ly = Ly_ef,
  Lt = Lt_ef,
  optns = optns
)

saveRDS(
  fpca_fit_ef,
  "path/fpca_fit_ef.rds"
)

# Use this instead when loading an existing result:
# fpca_fit_ef <- readRDS(
#   "path/fpca_fit_ef.rds"
# )

plot(fpca_fit_ef)

# Check that at least three FPCs were retained
stopifnot(ncol(fpca_fit_ef$xiEst) >= 3)

# ============================================================
# 7. Extract FPCA scores
# ============================================================

scores_ef <- as.data.frame(fpca_fit_ef$xiEst)

colnames(scores_ef) <- paste0(
  "FPC_EF",
  seq_len(ncol(scores_ef))
)

scores_ef$SubjID <- ids_ef

fpca_merged4 <- fpca_merged3 %>%
  select(-any_of(c(
    "baseline_EF",
    names(scores_ef)[names(scores_ef) != "SubjID"]
  ))) %>%
  left_join(scores_ef, by = "SubjID") %>%
  left_join(baseline_ef, by = "SubjID") %>%
  mutate(
    FPC_EF1_z = as.numeric(scale(FPC_EF1)),
    FPC_EF2_z = as.numeric(scale(FPC_EF2)),
  )

# ============================================================
# 8. Models predicting EF FPCA scores
# ============================================================

model_ef_fpc1 <- lm(
  FPC_EF1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total + phenoconversion_group,
  data = fpca_merged4
)

summary(model_ef_fpc1)

model_ef_fpc2 <- lm(
  FPC_EF2_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total+ phenoconversion_group,
  data = fpca_merged4
)

summary(model_ef_fpc2)

emtrends(
  model_ef_fpc1,
  ~ phenoconversion_group,
  var = "WMH_total"
)

# ============================================================
# 9. Plot EF FPCA components
# ============================================================

plot_fpc_effect_ef <- function(k) {

  lambda_k <- fpca_fit_ef$lambda[k]
  phi_k <- fpca_fit_ef$phi[, k]

  plot_df <- data.frame(
    age = fpca_fit_ef$workGrid,
    mean = fpca_fit_ef$mu,
    plus = fpca_fit_ef$mu +
      2 * sqrt(lambda_k) * phi_k,
    minus = fpca_fit_ef$mu -
      2 * sqrt(lambda_k) * phi_k
  )

  ggplot(plot_df, aes(x = age)) +
    geom_line(
      aes(y = mean),
      linewidth = 1.2
    ) +
    geom_line(
      aes(y = plus),
      linetype = "dotted",
      linewidth = 1
    ) +
    geom_line(
      aes(y = minus),
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste0("Executive function FPC", k),
      x = "Age at visit",
      y = "Executive function composite"
    ) +
    theme_classic()
}

plot_fpc_effect_ef(1)
plot_fpc_effect_ef(2)

# ============================================================
# 10. EF FPC scores by phenoconversion
# ============================================================

ggplot(
  fpca_merged4,
  aes(
    x = phenoconversion_group,
    y = FPC_EF1_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "Executive function FPC1 by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized EF FPC1 score"
  ) +
  theme_classic()

ggplot(
  fpca_merged4,
  aes(
    x = phenoconversion_group,
    y = FPC_EF2_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "Executive function FPC2 by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized EF FPC2 score"
  ) +
  theme_classic()

# ============================================================
# 11. Relationship with UPDRS FPCA
# ============================================================

ggplot(
  fpca_merged_ef,
  aes(
    x = FPC1,
    y = FPC_EF1_z,
    color = phenoconversion_group
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  labs(
    title = "UPDRS FPC1 versus executive-function FPC1",
    x = "UPDRS FPC1 score",
    y = "Standardized EF FPC1 score",
    color = "Phenoconversion status"
  ) +
  theme_classic()

model_ef_fpcx1 <- lm(
  FPC1 ~
    FPC_EF1_z +
    FPC_EF2_z +
    baseline_age +
    Sx +
    baseline_EF +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group,
  data = fpca_merged_ef
)

summary(model_ef_fpcx1)

model_ef_fpcx2 <- lm(
  FPC2 ~
    FPC_EF1_z +
    FPC_EF2_z +
    baseline_age +
    Sx +
    baseline_EF +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group,
  data = fpca_merged_ef
)

summary(model_ef_fpcx2)

# ============================================================
# 12. Logistic regression for phenoconversion
# ============================================================

fpca_merged_ef <- fpca_merged_ef %>%
  mutate(
    phenoconversion_group = factor(
      phenoconversion_group,
      levels = c(
        "Non-converter",
        "Phenoconverter"
      )
    )
  )

model_glm_ef <- glm(
  phenoconversion_group ~
    FPC1_z +
    FPC2_z +
    FPC_scopa1_z +
    FPC_scopa2_z +
    FPC_EF1_z +
    FPC_EF2_z +
    baseline_age +
    Sx +
    WMH_total,
  family = binomial(link = "logit"),
  data = fpca_merged4
)

summary(model_glm_ef)

########### Cognition M

dat_m <- prodromal_data %>%
  select(SubjID, Age, Cognition_M) %>%
  mutate(
    age_at_visit = as.numeric(Age),
    Memory = as.numeric(as.character(Cognition_M))
  ) %>%
  filter(
    !is.na(SubjID),
    !is.na(age_at_visit),
    !is.na(Memory)
  ) %>%
  arrange(SubjID, age_at_visit)

# ============================================================
# 2. Resolve duplicate visits at the same age
# ============================================================

dat_m <- dat_m %>%
  group_by(SubjID, age_at_visit) %>%
  summarise(
    Memory = mean(Memory, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SubjID, age_at_visit)

# Earliest available memory value for each participant
baseline_m <- dat_m %>%
  group_by(SubjID) %>%
  slice_min(
    order_by = age_at_visit,
    n = 1,
    with_ties = FALSE
  ) %>%
  transmute(
    SubjID,
    baseline_M = Memory
  ) %>%
  ungroup()

# ============================================================
# 3. Restrict age range
# ============================================================

dat2_m <- dat_m %>%
  filter(
    age_at_visit >= 50,
    age_at_visit <= 85
  )

# ============================================================
# 4. Create FPCA lists
# ============================================================

fpca_lists_m <- dat2_m %>%
  arrange(SubjID, age_at_visit) %>%
  group_by(SubjID) %>%
  summarise(
    Lt = list(age_at_visit),
    Ly = list(Memory),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  filter(n_visits >= 2)

ids_m <- fpca_lists_m$SubjID
Lt_m  <- fpca_lists_m$Lt
Ly_m  <- fpca_lists_m$Ly

# ============================================================
# 5. Check inputs
# ============================================================

stopifnot(length(ids_m) == length(Lt_m))
stopifnot(length(ids_m) == length(Ly_m))
stopifnot(all(lengths(Lt_m) == lengths(Ly_m)))

summary(lengths(Lt_m))
range(unlist(Lt_m), na.rm = TRUE)
range(unlist(Ly_m), na.rm = TRUE)

# ============================================================
# 6. Run FPCA
# ============================================================

optns <- list(
  dataType = "Sparse",
  kernel = "gauss",
  error = TRUE,
  methodSelectK = "FVE",
  FVEthreshold = 0.95,
  maxK = 5,
  verbose = TRUE
)

fpca_fit_m <- FPCA(
  Ly = Ly_m,
  Lt = Lt_m,
  optns = optns
)

saveRDS(
  fpca_fit_m,
  paste0(
    "path/fpca_fit_memory.rds"
  )
)

# Use this instead when loading an existing result:
# fpca_fit_m <- readRDS(
#   paste0(
#     "path/fpca_fit_memory.rds"
#   )
# )

plot(fpca_fit_m)

# Only two components are used below
stopifnot(ncol(fpca_fit_m$xiEst) >= 2)

# ============================================================
# 7. Extract FPCA scores
# ============================================================

scores_m <- as.data.frame(fpca_fit_m$xiEst)

colnames(scores_m) <- paste0(
  "FPC_M",
  seq_len(ncol(scores_m))
)

scores_m$SubjID <- ids_m

fpca_merged5 <- fpca_merged4 %>%
  select(
    -any_of(c(
      "baseline_M",
      "FPC_M1_z",
      "FPC_M2_z",
      names(scores_m)[names(scores_m) != "SubjID"]
    ))
  ) %>%
  left_join(scores_m, by = "SubjID") %>%
  left_join(baseline_m, by = "SubjID") %>%
  mutate(
    FPC_M1_z = as.numeric(scale(FPC_M1)),
    FPC_M2_z = as.numeric(scale(FPC_M2))
  )

# ============================================================
# 8. Models predicting memory FPCA scores
# ============================================================

model_m_fpc1 <- lm(
  FPC_M1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_m_fpc1)

model_m_fpc2 <- lm(
  FPC_M2_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_m_fpc2)

# ============================================================
# 8B. Test whether the WMH relationship differs by group
# ============================================================

model_m_fpc1_interaction <- lm(
  FPC_M1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total * phenoconversion_group,
  data = fpca_merged5
)

summary(model_m_fpc1_interaction)

# WMH slope within each phenoconversion group
emtrends(
  model_m_fpc1_interaction,
  ~ phenoconversion_group,
  var = "WMH_total"
)

# Compare the interaction model with the main-effects model
anova(
  model_m_fpc1,
  model_m_fpc1_interaction
)

# ============================================================
# 9. Plot memory FPCA components
# ============================================================

plot_fpc_effect_m <- function(k) {

  lambda_k <- fpca_fit_m$lambda[k]
  phi_k <- fpca_fit_m$phi[, k]

  plot_df <- data.frame(
    age = fpca_fit_m$workGrid,
    mean = fpca_fit_m$mu,
    plus = fpca_fit_m$mu +
      2 * sqrt(lambda_k) * phi_k,
    minus = fpca_fit_m$mu -
      2 * sqrt(lambda_k) * phi_k
  )

  ggplot(plot_df, aes(x = age)) +
    geom_line(
      aes(y = mean),
      linewidth = 1.2
    ) +
    geom_line(
      aes(y = plus),
      linetype = "dotted",
      linewidth = 1
    ) +
    geom_line(
      aes(y = minus),
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste0("Memory FPC", k),
      x = "Age at visit",
      y = "Memory composite"
    ) +
    theme_classic()
}

plot_fpc_effect_m(1)
plot_fpc_effect_m(2)

# ============================================================
# 10. Memory FPC scores by phenoconversion
# ============================================================

ggplot(
  fpca_merged5,
  aes(
    x = phenoconversion_group,
    y = FPC_M1_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "Memory FPC1 by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized memory FPC1 score"
  ) +
  theme_classic()

ggplot(
  fpca_merged5,
  aes(
    x = phenoconversion_group,
    y = FPC_M2_z
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "Memory FPC2 by phenoconversion status",
    x = "Phenoconversion status",
    y = "Standardized memory FPC2 score"
  ) +
  theme_classic()

# ============================================================
# 11. Relationship with UPDRS FPCA
# ============================================================

ggplot(
  fpca_merged5,
  aes(
    x = FPC1_z,
    y = FPC_M1_z,
    color = phenoconversion_group
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  labs(
    title = "UPDRS FPC1 versus memory FPC1",
    x = "Standardized UPDRS FPC1 score",
    y = "Standardized memory FPC1 score",
    color = "Phenoconversion status"
  ) +
  theme_classic()

model_m_fpcx1 <- lm(
  FPC1_z ~
    FPC_M1_z +
    FPC_M2_z +
    baseline_age +
    Sx +
    baseline_M +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_m_fpcx1)

model_m_fpcx2 <- lm(
  FPC2_z ~
    FPC_M1_z +
    FPC_M2_z +
    baseline_age +
    Sx +
    baseline_M +
    RBD +
    LRRK2 +
    GBA +
    hyposmia +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_m_fpcx2)

# ============================================================
# 12. Logistic regression for phenoconversion
# ============================================================

fpca_merged5 <- fpca_merged5 %>%
  mutate(
    phenoconversion_group = factor(
      phenoconversion_group,
      levels = c(
        "Non-converter",
        "Phenoconverter"
      )
    )
  )

model_glm_m <- glm(
  phenoconversion_group ~
    FPC1_z +
    FPC2_z +
    FPC_scopa1_z +
    FPC_scopa2_z +
    FPC_EF1_z +
    FPC_EF2_z +
    FPC_M1_z +
    FPC_M2_z +
    baseline_age +
    Sx + hyposmia + baseline_scopa + baseline_moca + baseline_EF + baseline_M +
    WMH_total,
  family = binomial(link = "logit"),
  data = fpca_merged5
)

summary(model_glm_m)

table(
  model.frame(model_glm_m)$phenoconversion_group
)

getwd()
write.csv(fpca_merged5,"path/FPCA_data.csv")

fpca_merged5 <- read.csv("path/FPCA_data.csv")

fpca_merged5 <- fpca_merged5 %>%
  mutate(
    phenoconversion_group = factor(
      phenoconversion_group,
      levels = c(
        "Non-converter",
        "Phenoconverter"
      )
    )
  )

fpca_merged5$WMH_total <- log1p(fpca_merged5$WMH_total)

########### All LMs here



model_fpc1 <- lm(
  FPC1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_fpc1)

model_fpc2 <- lm(
  FPC2_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_fpc2)

model_scopa_fpc1 <- lm(
  FPC_scopa1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total,
  data = fpca_merged5
)

summary(model_scopa_fpc1)

model_scopa_fpc2 <- lm(
  FPC_scopa2_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total,
  data = fpca_merged5
)

summary(model_scopa_fpc2)


model_EF_fpc1 <- lm(
  FPC_EF1_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total,
  data = fpca_merged5
)

summary(model_EF_fpc1)


model_EF_fpc2 <- lm(
  FPC_EF2_z ~
    baseline_age +
    Sx +
    orthostasis +
    hyposmia +
    WMH_total +
    phenoconversion_group,
  data = fpca_merged5
)

summary(model_EF_fpc2)


### function for FDA
baseline_chars <- prodromal_model_data %>%
  select(SubjID, years_from_baseline, Age, Sx, scopa,MDS_UPDRS3_OFF_SCORE,MoCA,Cognition_EF, Cognition_M,phenoconversion_group, hyposmia, orthostasis, RBD, LRRK2, GBA, SNCA, PRKN, Hyposmia,WMH_total,dWMH,pWMH,scaling_factor,Native_v) %>%
  mutate(
    Age = as.numeric(Age),
    baseline_updrs3 = as.numeric(as.character(MDS_UPDRS3_OFF_SCORE)),
    baseline_scopa = as.numeric(as.character(scopa)),
    baseline_moca = as.numeric(as.character(MoCA)),
    baseline_EF = as.numeric(as.character(Cognition_EF)),
    baseline_M = as.numeric(as.character(Cognition_M))


  ) %>%
  filter(!is.na(SubjID), !is.na(Age)) %>%
  arrange(SubjID, Age) %>%
  group_by(SubjID) %>%
  slice(1) %>%
  ungroup() %>%
  rename(
    baseline_age = Age
  ) 



# ============================================================
# SETTINGS
# ============================================================



MAX_FOLLOWUP <- 10.5
MIN_VISITS <- 2

# Change these only if your actual column names differ
OUTCOMES <- list(

  UPDRS = list(
    outcome = "MDS_UPDRS3_OFF_SCORE",
    prefix = "FPC",
    baseline = "baseline_updrs3"
  ),

  SCOPA = list(
    outcome = "scopa",
    prefix = "FPC_scopa",
    baseline = "baseline_scopa"
  ),

  MoCA = list(
    outcome = "MoCA",
    prefix = "FPC_moca",
    baseline = "baseline_moca"
  ),

  EF = list(
    outcome = "Cognition_EF",
    prefix = "FPC_EF",
    baseline = "baseline_EF"
  ),

  Memory = list(
    outcome = "Cognition_M",
    prefix = "FPC_M",
    baseline = "baseline_M"
  )
)


# Confirm that all outcome columns exist
required_outcomes <- vapply(
  OUTCOMES,
  function(x) x$outcome,
  character(1)
)

missing_outcomes <- setdiff(
  required_outcomes,
  names(prodromal_model_data)
)

if (length(missing_outcomes) > 0) {
  stop(
    "Missing outcome columns: ",
    paste(missing_outcomes, collapse = ", ")
  )
}


# ============================================================
# FPCA FUNCTION
# ============================================================

run_sparse_fpca <- function(
    data,
    outcome_column,
    prefix,
    max_followup = 10.5,
    min_visits = 2
) {

  # Prepare longitudinal data
  dat <- data %>%
    transmute(
      SubjID = as.character(SubjID),
      years_from_baseline =
        as.numeric(years_from_baseline),
      outcome =
        suppressWarnings(
          as.numeric(as.character(.data[[outcome_column]]))
        )
    ) %>%
    filter(
      !is.na(SubjID),
      is.finite(years_from_baseline),
      is.finite(outcome),
      years_from_baseline >= 0,
      years_from_baseline < max_followup
    ) %>%
    group_by(SubjID, years_from_baseline) %>%
    summarise(
      outcome = mean(outcome, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(SubjID, years_from_baseline)

  # Create lists for FPCA
  fpca_lists <- dat %>%
    group_by(SubjID) %>%
    summarise(
      Lt = list(years_from_baseline),
      Ly = list(outcome),
      n_visits = n(),
      .groups = "drop"
    ) %>%
    filter(n_visits >= min_visits)

  if (nrow(fpca_lists) == 0) {
    stop(
      "No eligible participants for ",
      outcome_column
    )
  }

  ids <- fpca_lists$SubjID
  Lt <- fpca_lists$Lt
  Ly <- fpca_lists$Ly

  stopifnot(length(ids) == length(Lt))
  stopifnot(length(ids) == length(Ly))
  stopifnot(all(lengths(Lt) == lengths(Ly)))

  # FPCA settings
  optns <- list(
    dataType = "Sparse",
    kernel = "gauss",
    error = TRUE,
    methodSelectK = "FVE",
    FVEthreshold = 0.95,
    maxK = 5,
    verbose = TRUE
  )

  # Run FPCA
  fit <- FPCA(
    Ly = Ly,
    Lt = Lt,
    optns = optns
  )

  if (ncol(fit$xiEst) < 2) {
    stop(
      outcome_column,
      " produced fewer than two FPCs."
    )
  }

  # Extract first two scores
  score_data <- data.frame(
    SubjID = ids,
    stringsAsFactors = FALSE
  )

  score_data[[paste0(prefix, "1")]] <-
    fit$xiEst[, 1]

  score_data[[paste0(prefix, "2")]] <-
    fit$xiEst[, 2]

  # Standardized scores
  score_data[[paste0(prefix, "1_z")]] <-
    as.numeric(
      scale(score_data[[paste0(prefix, "1")]])
    )

  score_data[[paste0(prefix, "2_z")]] <-
    as.numeric(
      scale(score_data[[paste0(prefix, "2")]])
    )

  list(
    fit = fit,
    scores = score_data,
    longitudinal_data = dat,
    Lt = Lt,
    Ly = Ly,
    visit_counts = fpca_lists %>%
      select(SubjID, n_visits)
  )
}


# ============================================================
# RUN ALL FIVE FPCAs
# ============================================================

fpca_results <- lapply(
  OUTCOMES,
  function(specification) {

    run_sparse_fpca(
      data = prodromal_model_data,
      outcome_column = specification$outcome,
      prefix = specification$prefix,
      max_followup = MAX_FOLLOWUP,
      min_visits = MIN_VISITS
    )
  }
)


# Access individual fits like this:
fpca_results$UPDRS$fit
fpca_results$SCOPA$fit
fpca_results$MoCA$fit
fpca_results$EF$fit
fpca_results$Memory$fit



plot(fpca_results$UPDRS$fit)
plot(fpca_results$SCOPA$fit)
plot(fpca_results$MoCA$fit)
plot(fpca_results$EF$fit)
plot(fpca_results$Memory$fit)


score_tables <- lapply(
  fpca_results,
  function(x) x$scores
)


baseline_for_merge <- baseline_chars %>%
  mutate(SubjID = as.character(SubjID))

score_tables <- lapply(
  score_tables,
  function(x) {
    x %>%
      mutate(SubjID = as.character(SubjID))
  }
)

fpca_merged5 <- Reduce(
  function(x, y) {
    left_join(x, y, by = "SubjID")
  },
  score_tables,
  init = baseline_for_merge
)

fpca_merged5 <- fpca_merged5 %>%
  mutate(
    phenoconversion_binary = case_when(
      phenoconversion_group %in% c("Converter", "Phenoconverter", "1", 1) ~ 1,
      phenoconversion_group %in% c("Non-converter", "Nonconverter", "0", 0) ~ 0,
      TRUE ~ NA_real_
    )
  )

plot_fpc_effect <- function(
    fit,
    k,
    outcome_name
) {

  lambda_k <- fit$lambda[k]
  phi_k <- fit$phi[, k]

  plot_df <- data.frame(
    time = fit$workGrid,
    mean = fit$mu,
    plus =
      fit$mu +
      2 * sqrt(lambda_k) * phi_k,
    minus =
      fit$mu -
      2 * sqrt(lambda_k) * phi_k
  )

  ggplot(plot_df, aes(x = time)) +
    geom_line(
      aes(y = mean),
      linewidth = 1.2
    ) +
    geom_line(
      aes(y = plus),
      linetype = "dotted",
      linewidth = 1
    ) +
    geom_line(
      aes(y = minus),
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste0(
        outcome_name,
        ": interpretation of FPC",
        k
      ),
      x = "Years from baseline",
      y = outcome_name
    ) +
    theme_classic()
}

fpca_results$SCOPA$fit

plot_fpc_effect(
  fpca_results$UPDRS$fit,
  2,
  "MDS-UPDRS III"
)

plot_fpc_effect(
  fpca_results$SCOPA$fit,
  1,
  "SCOPA-AUT"
)

plot_fpc_effect(
  fpca_results$MoCA$fit,
  2,
  "MoCA"
)

plot_fpc_effect(
  fpca_results$EF$fit,
  2,
  "Executive function"
)

plot_fpc_effect(
  fpca_results$Memory$fit,
  2,
  "Memory"
)


### FPC models with WMH
fit_fpc_model_family <- function(
    data,
    fpc1,
    fpc2,
    baseline_variable
) {

  adjustment_variables <- c(
    "baseline_age",
    "Sx",
    baseline_variable,
    "hyposmia",
    "orthostasis",
    "WMH_total"
  )

  required_variables <- c(
    "phenoconversion_binary",
    fpc1,
    fpc2,
    adjustment_variables
  )

  missing_variables <- setdiff(
    required_variables,
    names(data)
  )

  if (length(missing_variables) > 0) {
    stop(
      "Missing model variables: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  fpc1_formula <- reformulate(
    adjustment_variables,
    response = fpc1
  )

  fpc2_formula <- reformulate(
    adjustment_variables,
    response = fpc2
  )

  phenoconversion_formula <- reformulate(
    c(
      fpc1,
      fpc2,
      adjustment_variables
    ),
    response = "phenoconversion_binary"
  )

  list(

    FPC1_model = lm(
      fpc1_formula,
      data = data
    ),

    FPC2_model = lm(
      fpc2_formula,
      data = data
    ),

    phenoconversion_model = glm(
      phenoconversion_formula,
      data = data,
      family = binomial
    )
  )
}


### run models

all_models <- list(

  UPDRS = fit_fpc_model_family(
    data = fpca_merged5,
    fpc1 = "FPC1_z",
    fpc2 = "FPC2_z",
    baseline_variable = "baseline_scopa"
  ),

  SCOPA = fit_fpc_model_family(
    data = fpca_merged5,
    fpc1 = "FPC_scopa1_z",
    fpc2 = "FPC_scopa2_z",
    baseline_variable = "baseline_scopa"
  ),

  MoCA = fit_fpc_model_family(
    data = fpca_merged5,
    fpc1 = "FPC_moca1_z",
    fpc2 = "FPC_moca2_z",
    baseline_variable = "baseline_scopa"
  ),

  EF = fit_fpc_model_family(
    data = fpca_merged5,
    fpc1 = "FPC_EF1_z",
    fpc2 = "FPC_EF2_z",
    baseline_variable = "baseline_scopa"
  ),

  Memory = fit_fpc_model_family(
    data = fpca_merged5,
    fpc1 = "FPC_M1_z",
    fpc2 = "FPC_M2_z",
    baseline_variable = "baseline_scopa"
  )
)

model_summaries <- lapply(
  all_models,
  function(domain_models) {
    lapply(domain_models, summary)
  }
)

model_summaries



