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
prodromal_data <- "PATH OF FINAL SHEET" # Adjust to add sheet path


prodromal_data <- read.csv(prodromal_data)

## Descriptive details

prodromal_data <- prodromal_data %>%
  filter(Dx == 4)



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

mem_model <- lmer(
  MDS_UPDRS3_OFF_SCORE ~ baseline_age + Sx + Education_Years + years_from_baseline*phenoconversion_group*pWMH + Native_v + Equivalent_Vascular_RF +
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

summary(mem_model)

coef(summary(updrs_model))

confint(mem_model, method = "Wald")


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

# FPC functions with Plotting


# ============================================================
# BASELINE CHARACTERISTICS
# ============================================================

baseline_chars <- prodromal_model_data %>%
  select(
    SubjID,
    years_from_baseline,
    Age,
    Sx,
    scopa,
    MDS_UPDRS3_OFF_SCORE,
    MoCA,
    Cognition_EF,
    Cognition_M,
    phenoconversion_group,
    hyposmia,
    orthostasis,
    Education_Years,
    Equivalent_Vascular_RF,
    WMH_total,
    dWMH,
    pWMH,
    Native_v
  ) %>%
  mutate(
    Age = as.numeric(Age),
    baseline_updrs3 =
      as.numeric(as.character(MDS_UPDRS3_OFF_SCORE)),
    baseline_scopa =
      as.numeric(as.character(scopa)),
    baseline_moca =
      as.numeric(as.character(MoCA)),
    baseline_EF =
      as.numeric(as.character(Cognition_EF)),
    baseline_M =
      as.numeric(as.character(Cognition_M))
  ) %>%
  filter(
    !is.na(SubjID),
    is.finite(Age)
  ) %>%
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

MIN_AGE <- 50
MAX_AGE <- 85
MIN_VISITS <- 1


OUTCOMES <- list(

  UPDRS = list(
    outcome = "MDS_UPDRS3_OFF_SCORE",
    prefix = "FPC",
    label = "MDS-UPDRS III"
  ),

  SCOPA = list(
    outcome = "scopa",
    prefix = "FPC_scopa",
    label = "SCOPA-AUT"
  ),

  MoCA = list(
    outcome = "MoCA",
    prefix = "FPC_moca",
    label = "MoCA"
  ),

  EF = list(
    outcome = "Cognition_EF",
    prefix = "FPC_EF",
    label = "Executive function"
  ),

  Memory = list(
    outcome = "Cognition_M",
    prefix = "FPC_M",
    label = "Memory"
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
# AGE-BASED SPARSE FPCA FUNCTION
# ============================================================

run_sparse_fpca_age <- function(
    data,
    outcome_column,
    prefix,
    min_age = 50,
    max_age = 85,
    min_visits = 1
) {

  # Prepare longitudinal data using chronological age
  dat <- data %>%
    transmute(
      SubjID = as.character(SubjID),
      Age = suppressWarnings(
        as.numeric(as.character(Age))
      ),
      outcome = suppressWarnings(
        as.numeric(as.character(.data[[outcome_column]]))
      )
    ) %>%
    filter(
      !is.na(SubjID),
      is.finite(Age),
      is.finite(outcome),
      Age >= min_age,
      Age <= max_age
    ) %>%
    group_by(SubjID, Age) %>%
    summarise(
      outcome = mean(outcome, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(SubjID, Age)

  # Create subject-specific age and outcome lists
  fpca_lists <- dat %>%
    group_by(SubjID) %>%
    summarise(
      Lt = list(Age),
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

  optns <- list(
    dataType = "Sparse",
    kernel = "gauss",
    error = TRUE,
    methodSelectK = "FVE",
    FVEthreshold = 0.95,
    maxK = 5,
    verbose = TRUE
  )

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

  # Extract first two FPC scores
  score_data <- data.frame(
    SubjID = ids,
    stringsAsFactors = FALSE
  )

  score_data[[paste0(prefix, "1")]] <-
    fit$xiEst[, 1]

  score_data[[paste0(prefix, "2")]] <-
    fit$xiEst[, 2]

  # Standardize scores
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
# RUN ALL FIVE AGE-BASED FPCAs
# ============================================================

fpca_results <- lapply(
  OUTCOMES,
  function(specification) {

    run_sparse_fpca_age(
      data = prodromal_model_data,
      outcome_column = specification$outcome,
      prefix = specification$prefix,
      min_age = MIN_AGE,
      max_age = MAX_AGE,
      min_visits = MIN_VISITS
    )
  }
)


# Inspect individual fits
fpca_results$UPDRS$fit
fpca_results$SCOPA$fit
fpca_results$MoCA$fit
fpca_results$EF$fit
fpca_results$Memory$fit


# ============================================================
# MERGE FPC SCORES WITH BASELINE DATA
# ============================================================

score_tables <- lapply(
  fpca_results,
  function(x) x$scores
)

baseline_for_merge <- baseline_chars %>%
  mutate(
    SubjID = as.character(SubjID)
  )

score_tables <- lapply(
  score_tables,
  function(x) {
    x %>%
      mutate(
        SubjID = as.character(SubjID)
      )
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
      phenoconversion_group %in%
        c("Converter", "Phenoconverter", "1", 1) ~ 1,
      phenoconversion_group %in%
        c("Non-converter", "Nonconverter", "0", 0) ~ 0,
      TRUE ~ NA_real_
    )
  )


# ============================================================
# DRAFT INTERPRETATION OF EACH COMPONENT
# ============================================================
#
# This describes the direction of the reconstructed trajectory
# associated with a higher FPC score. The wording should still
# be checked against the plots before placing it in the thesis.
# ============================================================

draft_fpc_interpretation <- function(fit, k) {

  effect <- 1 * sqrt(fit$lambda[k]) * fit$phi[, k]

  tolerance <- 0.05 * max(abs(effect), na.rm = TRUE)

  n_grid <- length(effect)
  n_segment <- max(2, floor(0.20 * n_grid))

  younger_effect <- mean(
    effect[seq_len(n_segment)],
    na.rm = TRUE
  )

  older_effect <- mean(
    tail(effect, n_segment),
    na.rm = TRUE
  )

  direction <- function(x) {
    if (x > tolerance) {
      return("higher")
    }

    if (x < -tolerance) {
      return("lower")
    }

    "similar"
  }

  younger_direction <- direction(younger_effect)
  older_direction <- direction(older_effect)

  effect_sign <- ifelse(
    effect > tolerance,
    1,
    ifelse(effect < -tolerance, -1, 0)
  )

  effect_sign <- effect_sign[effect_sign != 0]

  number_crossings <- if (length(effect_sign) > 1) {
    sum(diff(effect_sign) != 0)
  } else {
    0
  }

  if (number_crossings > 1) {
    return(
      paste(
        "Complex age-varying pattern;",
        "interpretation should be based on the plotted trajectories"
      )
    )
  }

  if (younger_direction == older_direction) {

    if (younger_direction == "similar") {
      return(
        "Minimal separation from the mean trajectory across age"
      )
    }

    return(
      paste0(
        tools::toTitleCase(younger_direction),
        " outcome values across the observed age range"
      )
    )
  }

  if (younger_direction == "similar") {
    return(
      paste0(
        "Little separation at younger ages, followed by ",
        older_direction,
        " outcome values at older ages"
      )
    )
  }

  if (older_direction == "similar") {
    return(
      paste0(
        tools::toTitleCase(younger_direction),
        " outcome values at younger ages, with diminishing ",
        "separation at older ages"
      )
    )
  }

  paste0(
    tools::toTitleCase(younger_direction),
    " outcome values at younger ages and ",
    older_direction,
    " outcome values at older ages"
  )
}


# ============================================================
# FPCA SUMMARY TABLE
# ============================================================

extract_fpca_summary <- function(
    result,
    outcome_label
) {

  fit <- result$fit

  # Use cumulative FVE when available
  if (
    !is.null(fit$cumFVE) &&
    length(fit$cumFVE) >= 2
  ) {

    cumulative_fve <- as.numeric(fit$cumFVE)
    component_fve <- diff(c(0, cumulative_fve))

  } else {

    component_fve <-
      fit$lambda / sum(fit$lambda)

    cumulative_fve <- cumsum(component_fve)
  }

  tibble(
    Clinical_outcome = outcome_label,
    Participants_included = nrow(result$scores),
    FPC = c("FPC1", "FPC2"),
    Variance_explained_percent =
      round(100 * component_fve[1:2], 1),
    Cumulative_variance_explained_percent =
      round(100 * cumulative_fve[1:2], 1),
    Interpretation_of_higher_FPC_score = c(
      draft_fpc_interpretation(fit, 1),
      draft_fpc_interpretation(fit, 2)
    )
  )
}


fpca_summary_table <- bind_rows(
  lapply(
    names(fpca_results),
    function(outcome_key) {

      extract_fpca_summary(
        result = fpca_results[[outcome_key]],
        outcome_label = OUTCOMES[[outcome_key]]$label
      )
    }
  )
)

fpca_summary_table


# ============================================================
# CREATE PLOTTING DATA
# ============================================================

make_fpca_plot_data <- function(
    result,
    outcome_label
) {

  fit <- result$fit

  bind_rows(
    lapply(
      1:2,
      function(k) {

        effect_k <-
          1 * sqrt(fit$lambda[k]) * fit$phi[, k]

        tibble(
          Age = fit$workGrid,
          Outcome = outcome_label,
          FPC = paste0("FPC", k),
          Mean = fit$mu,
          `+1 SD` = fit$mu + effect_k,
          `-1 SD` = fit$mu - effect_k
        ) %>%
          pivot_longer(
            cols = c("Mean", "+1 SD", "-1 SD"),
            names_to = "Trajectory",
            values_to = "Outcome_value"
          )
      }
    )
  )
}


fpca_plot_data <- bind_rows(
  lapply(
    names(fpca_results),
    function(outcome_key) {

      make_fpca_plot_data(
        result = fpca_results[[outcome_key]],
        outcome_label = OUTCOMES[[outcome_key]]$label
      )
    }
  )
)


# Set panel order
fpca_plot_data <- fpca_plot_data %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(
        "MDS-UPDRS III",
        "SCOPA-AUT",
        "MoCA",
        "Executive function",
        "Memory"
      )
    ),
    FPC = factor(
      FPC,
      levels = c("FPC1", "FPC2")
    ),
    Trajectory = factor(
      Trajectory,
      levels = c("Mean", "+1 SD", "-1 SD")
    )
  )


# ============================================================
# COMBINED 5 × 2 FPCA FIGURE
# ============================================================

fpca_combined_plot <- ggplot(
  fpca_plot_data,
  aes(
    x = Age,
    y = Outcome_value,
    colour = Trajectory
  )
) +
  geom_line(
    linewidth = 1
  ) +
  facet_grid(
    rows = vars(Outcome),
    cols = vars(FPC),
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "Mean" = "black",
      "+1 SD" = "#0072B2",
      "-1 SD" = "#D62728"
    ),
    labels = c(
      "Mean trajectory",
      "+1 SD (higher FPC score)",
      "-1 SD (lower FPC score)"
    )
  ) +
  labs(
    x = "Age, years",
    y = "Estimated outcome value",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey70"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(1, "lines")
  )

fpca_combined_plot


# Optional high-resolution export
ggsave(
  filename = "FPCA_age_based_trajectories.png",
  plot = fpca_combined_plot,
  width = 10,
  height = 14,
  units = "in",
  dpi = 600
)


# ============================================================
# SEPARATE TWO-PANEL PLOTS
# ============================================================

plot_fpca_outcome <- function(outcome_name) {

  fpca_plot_data %>%
    filter(Outcome == outcome_name) %>%
    ggplot(
      aes(
        x = Age,
        y = Outcome_value,
        colour = Trajectory
      )
    ) +
    geom_line(
      linewidth = 1.1
    ) +
    facet_wrap(
      ~FPC,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_colour_manual(
      values = c(
        "Mean" = "black",
        "+1 SD" = "#0072B2",
        "-1 SD" = "#D62728"
      ),
      labels = c(
        "Mean trajectory",
        "+1 SD (higher FPC score)",
        "-1 SD (lower FPC score)"
      )
    ) +
    labs(
      title = outcome_name,
      x = "Age, years",
      y = "Estimated outcome value",
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey70"
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
}


updrs_fpca_plot <- plot_fpca_outcome("MDS-UPDRS III")
scopa_fpca_plot <- plot_fpca_outcome("SCOPA-AUT")
moca_fpca_plot <- plot_fpca_outcome("MoCA")
ef_fpca_plot <- plot_fpca_outcome("Executive function")
memory_fpca_plot <- plot_fpca_outcome("Memory")


updrs_fpca_plot
scopa_fpca_plot
moca_fpca_plot
ef_fpca_plot
memory_fpca_plot



# ============================================================
# MODEL SETTINGS
# ============================================================

# Demographic adjustments selected a priori
DEMOGRAPHIC_COVARIATES <- c(
  "baseline_age",
  "Sx"
)

# Additional clinical covariates
# Remove either variable if it is not part of your prespecified model
ADDITIONAL_COVARIATES <- c(
  "Education_Years",
  "Equivalent_Vascular_RF",
  "WMH_total",
  "Native_v"
)

COMMON_ADJUSTMENTS <- c(
  DEMOGRAPHIC_COVARIATES,
  ADDITIONAL_COVARIATES
)

# Imaging exposure of interest
WMH_EXPOSURE <- "WMH_total"


# ============================================================
# PREPARE BASELINE MODEL DATA
# ============================================================

baseline_model_data <- baseline_chars %>%
  mutate(
    SubjID = as.character(SubjID),

    phenoconversion_character = tolower(
      trimws(as.character(phenoconversion_group))
    ),

    phenoconversion_binary = case_when(
      phenoconversion_character %in%
        c("converter", "phenoconverter", "1") ~ 1,

      phenoconversion_character %in%
        c("non-converter", "nonconverter", "0") ~ 0,

      TRUE ~ NA_real_
    )
  ) %>%
  select(-phenoconversion_character)


# ============================================================
# FUNCTION FOR ONE OUTCOME
# ============================================================

fit_fpc_model_family <- function(
    result,
    baseline_data,
    fpc1,
    fpc2,
    baseline_variable,
    common_adjustments,
    exposure = "WMH_total",
    min_visits_for_models = 1
) {

  # Join the domain-specific scores and visit counts
  model_data <- baseline_data %>%
    inner_join(
      result$scores,
      by = "SubjID"
    ) %>%
    inner_join(
      result$visit_counts,
      by = "SubjID"
    ) %>%
    filter(
      n_visits >= min_visits_for_models
    )

  required_variables <- unique(
    c(
      "phenoconversion_binary",
      "n_visits",
      fpc1,
      fpc2,
      exposure,
      baseline_variable,
      common_adjustments
    )
  )

  missing_variables <- setdiff(
    required_variables,
    names(model_data)
  )

  if (length(missing_variables) > 0) {
    stop(
      "Missing model variables: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  # Models of clinical trajectory
  #
  # Each model includes:
  #   1. WMH burden
  #   2. Baseline value of the corresponding clinical outcome
  #   3. Selected demographic and clinical covariates

  fpc_predictors <- unique(
    c(
      exposure,
      baseline_variable,
      common_adjustments
    )
  )

  fpc1_formula <- reformulate(
    termlabels = fpc_predictors,
    response = fpc1
  )

  fpc2_formula <- reformulate(
    termlabels = fpc_predictors,
    response = fpc2
  )

  # Phenoconversion model
  #
  # Tests whether the FPCA-derived trajectories are associated
  # with phenoconversion after adjustment for WMH, baseline
  # clinical severity, and the selected covariates.

  phenoconversion_predictors <- unique(
    c(
      fpc1,
      fpc2,
      exposure,
      baseline_variable,
      common_adjustments
    )
  )

  phenoconversion_formula <- reformulate(
    termlabels = phenoconversion_predictors,
    response = "phenoconversion_binary"
  )

  fpc1_model <- lm(
    formula = fpc1_formula,
    data = model_data,
    na.action = na.exclude
  )

  fpc2_model <- lm(
    formula = fpc2_formula,
    data = model_data,
    na.action = na.exclude
  )

  phenoconversion_model <- glm(
    formula = phenoconversion_formula,
    data = model_data,
    family = binomial,
    na.action = na.exclude
  )

  list(
    model_data = model_data,

    eligible_participants = nrow(model_data),

    FPC1_model = fpc1_model,
    FPC2_model = fpc2_model,
    phenoconversion_model = phenoconversion_model,

    model_sample_sizes = c(
      FPC1 = nobs(fpc1_model),
      FPC2 = nobs(fpc2_model),
      Phenoconversion = nobs(phenoconversion_model)
    )
  )
}


# ============================================================
# DOMAIN-SPECIFIC MODEL SPECIFICATIONS
# ============================================================

MODEL_SPECIFICATIONS <- list(

  UPDRS = list(
    fpc1 = "FPC1_z",
    fpc2 = "FPC2_z",
    baseline = "baseline_updrs3"
  ),

  SCOPA = list(
    fpc1 = "FPC_scopa1_z",
    fpc2 = "FPC_scopa2_z",
    baseline = "baseline_scopa"
  ),

  MoCA = list(
    fpc1 = "FPC_moca1_z",
    fpc2 = "FPC_moca2_z",
    baseline = "baseline_moca"
  ),

  EF = list(
    fpc1 = "FPC_EF1_z",
    fpc2 = "FPC_EF2_z",
    baseline = "baseline_EF"
  ),

  Memory = list(
    fpc1 = "FPC_M1_z",
    fpc2 = "FPC_M2_z",
    baseline = "baseline_M"
  )
)


# ============================================================
# FUNCTION TO RUN ALL MODEL FAMILIES
# ============================================================

run_all_fpc_models <- function(
    minimum_visits
) {

  setNames(
    lapply(
      names(MODEL_SPECIFICATIONS),
      function(domain) {

        specification <-
          MODEL_SPECIFICATIONS[[domain]]

        fit_fpc_model_family(
          result = fpca_results[[domain]],
          baseline_data = baseline_model_data,

          fpc1 = specification$fpc1,
          fpc2 = specification$fpc2,
          baseline_variable = specification$baseline,

          common_adjustments = COMMON_ADJUSTMENTS,
          exposure = WMH_EXPOSURE,

          min_visits_for_models = minimum_visits
        )
      }
    ),
    names(MODEL_SPECIFICATIONS)
  )
}


# ============================================================
# PRIMARY ANALYSIS: AT LEAST TWO ASSESSMENTS
# ============================================================

primary_models <- run_all_fpc_models(
  minimum_visits = 1
)



# ============================================================
# VIEW INDIVIDUAL MODEL RESULTS
# ============================================================

# MDS-UPDRS III
summary(primary_models$UPDRS$FPC1_model)
summary(primary_models$UPDRS$FPC2_model)
summary(primary_models$UPDRS$phenoconversion_model)

# SCOPA-AUT
summary(primary_models$SCOPA$FPC1_model)
summary(primary_models$SCOPA$FPC2_model)
summary(primary_models$SCOPA$phenoconversion_model)

# MoCA
summary(primary_models$MoCA$FPC1_model)
summary(primary_models$MoCA$FPC2_model)
summary(primary_models$MoCA$phenoconversion_model)

# Executive function
summary(primary_models$EF$FPC1_model)
summary(primary_models$EF$FPC2_model)
summary(primary_models$EF$phenoconversion_model)

# Memory
summary(primary_models$Memory$FPC1_model)
summary(primary_models$Memory$FPC2_model)
summary(primary_models$Memory$phenoconversion_model)


# ============================================================
# SAMPLE SIZE USED BY EACH MODEL
# ============================================================

primary_model_sample_sizes <- bind_rows(
  lapply(
    names(primary_models),
    function(domain) {

      tibble(
        Clinical_outcome = domain,
        Eligible_with_required_visits =
          primary_models[[domain]]$eligible_participants,

        FPC1_model_N =
          primary_models[[domain]]$model_sample_sizes["FPC1"],

        FPC2_model_N =
          primary_models[[domain]]$model_sample_sizes["FPC2"],

        Phenoconversion_model_N =
          primary_models[[domain]]$
            model_sample_sizes["Phenoconversion"]
      )
    }
  )
)

primary_model_sample_sizes




## Prepare table

# ============================================================
# TABLE: ASSOCIATION BETWEEN WMH_total AND FPC SCORES
# PRIMARY ANALYSIS: PARTICIPANTS WITH >= 2 VISITS
# ============================================================


# Labels used in the final table
OUTCOME_LABELS <- c(
  UPDRS = "MDS-UPDRS III",
  SCOPA = "SCOPA-AUT",
  MoCA = "MoCA",
  EF = "Executive function",
  Memory = "Memory"
)


# Extract WMH_total coefficient from one linear model
extract_wmh_result <- function(
    model,
    clinical_outcome,
    component,
    exposure = "WMH_total"
) {

  coefficient_table <- summary(model)$coefficients

  if (!exposure %in% rownames(coefficient_table)) {
    stop(
      exposure,
      " was not found in the ",
      clinical_outcome,
      " ",
      component,
      " model."
    )
  }

  beta <- coefficient_table[exposure, "Estimate"]
  standard_error <- coefficient_table[exposure, "Std. Error"]
  p_value <- coefficient_table[exposure, "Pr(>|t|)"]

  confidence_interval <- confint(
    model,
    parm = exposure,
    level = 0.95
  )

  tibble(
    Clinical_outcome = clinical_outcome,
    Component = component,
    N = nobs(model),
    Beta = beta,
    Standard_error = standard_error,
    CI_lower = confidence_interval[1],
    CI_upper = confidence_interval[2],
    P_value = p_value
  )
}


# Extract results from all primary models
primary_wmh_results <- map_dfr(
  names(primary_models),
  function(domain) {

    domain_models <- primary_models[[domain]]
    outcome_label <- OUTCOME_LABELS[[domain]]

    bind_rows(
      extract_wmh_result(
        model = domain_models$FPC1_model,
        clinical_outcome = outcome_label,
        component = "FPC1"
      ),

      extract_wmh_result(
        model = domain_models$FPC2_model,
        clinical_outcome = outcome_label,
        component = "FPC2"
      )
    )
  }
)


# ============================================================
# FORMAT AS A DISSERTATION-READY TABLE
# ============================================================

primary_wmh_table <- primary_wmh_results %>%
  mutate(
    `β (95% CI)` = sprintf(
      "%.3f (%.3f to %.3f)",
      Beta,
      CI_lower,
      CI_upper
    ),

    `P value` = case_when(
      P_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_value)
    )
  ) %>%
  select(
    `Clinical outcome` = Clinical_outcome,
    `Functional component` = Component,
    N,
    `β (95% CI)`,
    `P value`
  )


primary_wmh_table

primary_wmh_table_fdr <- primary_wmh_results %>%
  mutate(
    FDR_adjusted_P = p.adjust(
      P_value,
      method = "BH"
    ),

    `β (95% CI)` = sprintf(
      "%.3f (%.3f to %.3f)",
      Beta,
      CI_lower,
      CI_upper
    ),

    `P value` = ifelse(
      P_value < 0.001,
      "<0.001",
      sprintf("%.3f", P_value)
    ),

    `FDR-adjusted P value` = ifelse(
      FDR_adjusted_P < 0.001,
      "<0.001",
      sprintf("%.3f", FDR_adjusted_P)
    )
  ) %>%
  select(
    `Clinical outcome` = Clinical_outcome,
    `Functional component` = Component,
    N,
    `β (95% CI)`,
    `P value`,
    `FDR-adjusted P value`
  )

primary_wmh_table_fdr



### Make more tables

# ============================================================
# WMH EXPOSURES
# ============================================================

# These are fitted in separate models
WMH_EXPOSURES <- c(
  WMH_total = "Total WMH",
  pWMH = "Periventricular WMH",
  dWMH = "Deep WMH"
)


# ============================================================
# CLINICAL OUTCOME LABELS
# ============================================================

OUTCOME_LABELS <- c(
  UPDRS = "MDS-UPDRS III",
  SCOPA = "SCOPA-AUT",
  MoCA = "MoCA",
  EF = "Executive function",
  Memory = "Memory"
)


# ============================================================
# FPCA MODEL SPECIFICATIONS
# ============================================================

FPCA_MODEL_SPECIFICATIONS <- list(

  UPDRS = list(
    fpc1 = "FPC1_z",
    fpc2 = "FPC2_z",
    baseline = "baseline_updrs3"
  ),

  SCOPA = list(
    fpc1 = "FPC_scopa1_z",
    fpc2 = "FPC_scopa2_z",
    baseline = "baseline_scopa"
  ),

  MoCA = list(
    fpc1 = "FPC_moca1_z",
    fpc2 = "FPC_moca2_z",
    baseline = "baseline_moca"
  ),

  EF = list(
    fpc1 = "FPC_EF1_z",
    fpc2 = "FPC_EF2_z",
    baseline = "baseline_EF"
  ),

  Memory = list(
    fpc1 = "FPC_M1_z",
    fpc2 = "FPC_M2_z",
    baseline = "baseline_M"
  )
)


# Do not include WMH_total, pWMH or dWMH here.
# The exposure will be added separately by the model function.

FPCA_COMMON_ADJUSTMENTS <- c(
  "baseline_age",
  "Sx",
  "Education_Years",
  "Equivalent_Vascular_RF",
  "Native_v"
)

# ============================================================
# FIT ONE FPCA MODEL FAMILY
# ============================================================

fit_fpca_wmh_models <- function(
    result,
    baseline_data,
    fpc1,
    fpc2,
    baseline_variable,
    exposure,
    adjustments,
    minimum_visits = 1
) {

  model_data <- baseline_data %>%
    mutate(
      SubjID = as.character(SubjID)
    ) %>%
    inner_join(
      result$scores %>%
        mutate(SubjID = as.character(SubjID)),
      by = "SubjID"
    ) %>%
    inner_join(
      result$visit_counts %>%
        mutate(SubjID = as.character(SubjID)),
      by = "SubjID"
    ) %>%
    filter(
      n_visits >= minimum_visits
    )

  predictors <- unique(
    c(
      exposure,
      baseline_variable,
      adjustments
    )
  )

  required_variables <- unique(
    c(
      fpc1,
      fpc2,
      predictors
    )
  )

  missing_variables <- setdiff(
    required_variables,
    names(model_data)
  )

  if (length(missing_variables) > 0) {
    stop(
      "Missing variables: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  fpc1_model <- lm(
    reformulate(
      termlabels = predictors,
      response = fpc1
    ),
    data = model_data,
    na.action = na.exclude
  )

  fpc2_model <- lm(
    reformulate(
      termlabels = predictors,
      response = fpc2
    ),
    data = model_data,
    na.action = na.exclude
  )

  list(
    model_data = model_data,
    FPC1_model = fpc1_model,
    FPC2_model = fpc2_model
  )
}


# ============================================================
# RUN 5 DOMAINS × 3 WMH EXPOSURES
# ============================================================

primary_fpca_models <- imap(
  FPCA_MODEL_SPECIFICATIONS,
  function(specification, domain) {

    imap(
      WMH_EXPOSURES,
      function(exposure_label, exposure_variable) {

        fit_fpca_wmh_models(
          result = fpca_results[[domain]],
          baseline_data = baseline_model_data,

          fpc1 = specification$fpc1,
          fpc2 = specification$fpc2,
          baseline_variable = specification$baseline,

          exposure = exposure_variable,
          adjustments = FPCA_COMMON_ADJUSTMENTS,
          minimum_visits = 1
        )
      }
    )
  }
)

# ============================================================
# EXTRACT ONE FPCA EXPOSURE RESULT
# ============================================================

extract_fpca_exposure <- function(
    model,
    clinical_outcome,
    component,
    exposure,
    exposure_label
) {

  coefficient_table <- summary(model)$coefficients

  if (!exposure %in% rownames(coefficient_table)) {
    stop(
      exposure,
      " not found in the model."
    )
  }

  confidence_interval <- confint(
    model,
    parm = exposure,
    level = 0.95
  )

  tibble(
    Analysis = "FPCA",
    Clinical_outcome = clinical_outcome,
    Exposure = exposure_label,
    Exposure_variable = exposure,
    Result_type = component,
    Association = paste0(component, " score"),
    Participants = nobs(model),
    Observations = NA_integer_,
    Beta = coefficient_table[exposure, "Estimate"],
    Standard_error =
      coefficient_table[exposure, "Std. Error"],
    CI_lower = confidence_interval[1],
    CI_upper = confidence_interval[2],
    P_value =
      coefficient_table[exposure, "Pr(>|t|)"]
  )
}


# ============================================================
# CREATE COMPLETE FPCA RESULTS TABLE
# ============================================================

fpca_wmh_results <- imap_dfr(
  primary_fpca_models,
  function(domain_models, domain) {

    imap_dfr(
      domain_models,
      function(exposure_models, exposure_variable) {

        exposure_label <-
          unname(WMH_EXPOSURES[exposure_variable])

        bind_rows(

          extract_fpca_exposure(
            model = exposure_models$FPC1_model,
            clinical_outcome = OUTCOME_LABELS[[domain]],
            component = "FPC1",
            exposure = exposure_variable,
            exposure_label = exposure_label
          ),

          extract_fpca_exposure(
            model = exposure_models$FPC2_model,
            clinical_outcome = OUTCOME_LABELS[[domain]],
            component = "FPC2",
            exposure = exposure_variable,
            exposure_label = exposure_label
          )
        )
      }
    )
  }
)

fpca_wmh_results

# ============================================================
# PREPARE MIXED-EFFECTS DATA
# ============================================================

mixed_model_data <- prodromal_model_data %>%
  mutate(
    SubjID = factor(SubjID),

    phenoconversion_group = relevel(
      factor(phenoconversion_group),
      ref = "Non-converter"
    )
  )

levels(mixed_model_data$phenoconversion_group)


# ============================================================
# REPEATED CLINICAL OUTCOMES
# ============================================================

MIXED_MODEL_OUTCOMES <- c(
  UPDRS = "MDS_UPDRS3_OFF_SCORE",
  SCOPA = "scopa",
  MoCA = "MoCA",
  EF = "Cognition_EF",
  Memory = "Cognition_M"
)


# ============================================================
# FIT ONE MIXED-EFFECTS MODEL
# ============================================================

fit_one_mixed_wmh_model <- function(
    data,
    outcome,
    exposure
) {

  model_formula <- as.formula(
    paste0(
      outcome,
      " ~ ",
      "Age + ",
      "Sx + ",
      "Education_Years + ",
      "phenoconversion_group * ",
      exposure,
      " + ",
      "Native_v + ",
      "Equivalent_Vascular_RF + ",
      "(1 | SubjID)"
    )
  )

  lmer(
    formula = model_formula,
    data = data,
    REML = FALSE,
    na.action = na.exclude,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
}


# ============================================================
# RUN 5 OUTCOMES × 3 WMH EXPOSURES
# ============================================================

mixed_wmh_models <- imap(
  MIXED_MODEL_OUTCOMES,
  function(outcome_variable, domain) {

    imap(
      WMH_EXPOSURES,
      function(exposure_label, exposure_variable) {

        fit_one_mixed_wmh_model(
          data = mixed_model_data,
          outcome = outcome_variable,
          exposure = exposure_variable
        )
      }
    )
  }
)

###

# ============================================================
# EXTRACT MIXED-EFFECTS WMH RESULTS
# ============================================================

extract_mixed_wmh_results <- function(
    model,
    clinical_outcome,
    exposure,
    exposure_label,
    reference_group = "Non-converter",
    comparison_group = "Converter"
) {

  beta <- fixef(model)
  variance_covariance <- as.matrix(vcov(model))
  coefficient_names <- names(beta)

  if (!exposure %in% coefficient_names) {
    stop(
      "Main exposure term ",
      exposure,
      " not found."
    )
  }

  interaction_term <- coefficient_names[
    grepl(":", coefficient_names, fixed = TRUE) &
      grepl(exposure, coefficient_names, fixed = TRUE) &
      grepl(
        "phenoconversion_group",
        coefficient_names,
        fixed = TRUE
      )
  ]

  if (length(interaction_term) != 1) {
    stop(
      "Expected one interaction term for ",
      exposure,
      "; found: ",
      paste(interaction_term, collapse = ", ")
    )
  }

  model_frame <- model.frame(model)

  participants <- n_distinct(
    model_frame$SubjID
  )

  observations <- nobs(model)

  make_result_row <- function(
      result_type,
      association,
      estimate,
      variance
  ) {

    standard_error <- sqrt(variance)
    z_value <- estimate / standard_error

    p_value <- 2 * pnorm(
      abs(z_value),
      lower.tail = FALSE
    )

    tibble(
      Analysis = "Linear mixed-effects model",
      Clinical_outcome = clinical_outcome,
      Exposure = exposure_label,
      Exposure_variable = exposure,
      Result_type = result_type,
      Association = association,
      Participants = participants,
      Observations = observations,
      Beta = estimate,
      Standard_error = standard_error,
      CI_lower =
        estimate - 1.96 * standard_error,
      CI_upper =
        estimate + 1.96 * standard_error,
      P_value = p_value
    )
  }

  # Association in the reference group
  reference_result <- make_result_row(
    result_type = "Non-converter slope",

    association = paste0(
      "WMH association: ",
      reference_group
    ),

    estimate = beta[exposure],

    variance =
      variance_covariance[exposure, exposure]
  )

  # Formal difference between groups
  interaction_result <- make_result_row(
    result_type = "Interaction",

    association = paste0(
      "Difference in WMH association: ",
      comparison_group,
      " vs ",
      reference_group
    ),

    estimate = beta[interaction_term],

    variance =
      variance_covariance[
        interaction_term,
        interaction_term
      ]
  )

  # Converter slope = main effect + interaction
  comparison_estimate <-
    beta[exposure] +
    beta[interaction_term]

  comparison_variance <-
    variance_covariance[exposure, exposure] +
    variance_covariance[
      interaction_term,
      interaction_term
    ] +
    2 * variance_covariance[
      exposure,
      interaction_term
    ]

  comparison_result <- make_result_row(
    result_type = "Converter slope",

    association = paste0(
      "WMH association: ",
      comparison_group
    ),

    estimate = comparison_estimate,
    variance = comparison_variance
  )

  bind_rows(
    reference_result,
    interaction_result,
    comparison_result
  )
}


# ============================================================
# CREATE COMPLETE MIXED-EFFECTS RESULTS
# ============================================================

mixed_wmh_results <- imap_dfr(
  mixed_wmh_models,
  function(domain_models, domain) {

    imap_dfr(
      domain_models,
      function(model, exposure_variable) {

        extract_mixed_wmh_results(
          model = model,
          clinical_outcome =
            OUTCOME_LABELS[[domain]],

          exposure = exposure_variable,

          exposure_label =
            unname(
              WMH_EXPOSURES[exposure_variable]
            )
        )
      }
    )
  }
)

mixed_wmh_results


# ============================================================
# COMBINE FPCA AND MIXED-EFFECTS RESULTS
# ============================================================

comprehensive_wmh_results <- bind_rows(
  fpca_wmh_results,
  mixed_wmh_results
) %>%
  group_by(
    Analysis,
    Result_type
  ) %>%
  mutate(
    FDR_adjusted_P = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  ungroup()


# ============================================================
# FORMAT FOR REPORTING
# ============================================================

comprehensive_wmh_table <- comprehensive_wmh_results %>%
  mutate(
    Panel = case_when(
      Analysis == "FPCA" ~
        "A. Functional principal component analysis",

      Analysis == "Linear mixed-effects model" ~
        "B. Longitudinal mixed-effects analysis"
    ),

    Sample = case_when(
      is.na(Observations) ~
        paste0(
          Participants,
          " participants"
        ),

      TRUE ~
        paste0(
          Participants,
          " participants; ",
          Observations,
          " observations"
        )
    ),

    `β (95% CI)` = sprintf(
      "%.3f (%.3f to %.3f)",
      Beta,
      CI_lower,
      CI_upper
    ),

    `P value` = case_when(
      P_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_value)
    ),

    `FDR-adjusted P value` = case_when(
      FDR_adjusted_P < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        FDR_adjusted_P
      )
    )
  ) %>%
  select(
    Panel,
    `Clinical outcome` = Clinical_outcome,
    Exposure,
    Association,
    Sample,
    `β (95% CI)`,
    `P value`,
    `FDR-adjusted P value`
  ) %>%
  arrange(
    Panel,
    `Clinical outcome`,
    Exposure,
    Association
  )

comprehensive_wmh_table


write.csv(
  comprehensive_wmh_table,
  "Comprehensive_global_WMH_results_formatted.csv",
  row.names = FALSE
)

write.csv(
  comprehensive_wmh_results,
  "Comprehensive_global_WMH_results_numeric.csv",
  row.names = FALSE
)
