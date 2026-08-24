# ============================================================================
# GLOBAL WMH MODELS AFTER FPCA: ALL AVAILABLE ASSESSMENTS
#
# Run this script only AFTER your existing FPCA code has created:
#   1. prodromal_long
#   2. fpca_results
#
# This script DOES NOT rerun FPCA.
#
# It produces three model panels:
#   A. Existing FPC scores ~ baseline WMH * phenoconversion group
#   B. Repeated clinical outcomes ~ time * baseline WMH * group + covariates
#   C. Repeated clinical outcomes ~ time * visit-matched WMH * group
#      + covariates
#
# Important decisions:
#   - FPCA and baseline-WMH analyses do not impose a minimum number of
#     repeated assessments.
#   - Visit-matched longitudinal-WMH models require at least two complete,
#     matched WMH-clinical visits per participant in each fitted model.
#   - No within-person/between-person decomposition is performed.
#   - Baseline means EVENT_ID == "BL".
#   - Longitudinal WMH is matched to clinical data by SubjID + EVENT_ID.
#   - Total, periventricular and deep WMH are fitted in separate models.
#   - Model sample sizes come from the complete-case rows actually fitted.
# ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lme4)


# ============================================================================
# 1. SETTINGS
# ============================================================================

BASELINE_EVENT <- "BL"
MIN_AGE <- 50
MAX_AGE <- 85

# These source columns match the later modelling section of your supplied code.
# The analysis variables are log(1 + native ROI WMH volume).
WMH_SOURCE_COLUMNS <- c(
  WMH_total = "WMH_native_roi",
  pWMH = "WMH_peri_roi",
  dWMH = "WMH_deep_roi"
)

WMH_LABELS <- c(
  WMH_total = "Total WMH",
  pWMH = "Periventricular WMH",
  dWMH = "Deep WMH"
)

CLINICAL_OUTCOMES <- c(
  UPDRS = "MDS_UPDRS3_OFF_SCORE",
  SCOPA = "scopa",
  MoCA = "MoCA",
  EF = "Cognition_EF",
  Memory = "Cognition_M"
)

OUTCOME_LABELS <- c(
  UPDRS = "MDS-UPDRS III",
  SCOPA = "SCOPA-AUT",
  MoCA = "MoCA",
  EF = "Executive function",
  Memory = "Memory"
)

FPCA_SPECIFICATIONS <- list(
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

# Baseline covariates used in all models.
# Native brain volume is added separately because it is baseline-specific in
# panels A/B and visit-specific in panel C.
MODEL_ADJUSTMENTS <- c(
  "baseline_age",
  "Sx",
  "Education_Years",
  "Equivalent_Vascular_RF"
)


# ============================================================================
# 2. VALIDATION AND FORMATTING HELPERS
# ============================================================================

assert_columns <- function(data, variables, data_name) {
  missing_variables <- setdiff(variables, names(data))

  if (length(missing_variables) > 0) {
    stop(
      data_name,
      " is missing: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  invisible(data)
}

assert_unique_keys <- function(data, keys, data_name) {
  duplicated_keys <- data %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)

  if (nrow(duplicated_keys) > 0) {
    print(head(duplicated_keys, 20))
    stop(
      data_name,
      " contains duplicated ",
      paste(keys, collapse = " + "),
      " keys. The script stopped to prevent a many-to-many merge."
    )
  }

  invisible(data)
}

assert_no_row_expansion <- function(before_rows, after_rows, join_name) {
  if (after_rows > before_rows) {
    stop(
      join_name,
      " expanded rows from ",
      before_rows,
      " to ",
      after_rows,
      ". Check the join keys."
    )
  }
}

numeric_clean <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_log1p <- function(x) {
  x <- numeric_clean(x)
  ifelse(is.finite(x) & x >= 0, log1p(x), NA_real_)
}

normalise_group <- function(x) {
  x_character <- tolower(trimws(as.character(x)))

  case_when(
    x_character %in% c(
      "1", "true", "yes", "converter", "phenoconverter",
      "phenoconvertor"
    ) ~ "Phenoconverter",
    x_character %in% c(
      "0", "false", "no", "non-converter", "nonconverter",
      "non-convertor", "nonconvertor"
    ) ~ "Non-converter",
    TRUE ~ NA_character_
  )
}

complete_model_data <- function(data, variables) {
  data %>%
    select(all_of(unique(variables))) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
}

format_p <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", x)
  )
}

fixed_effects <- function(model) {
  if (inherits(model, "merMod")) {
    lme4::fixef(model)
  } else {
    coef(model)
  }
}

find_exact_term <- function(model, required_parts) {
  coefficient_names <- names(fixed_effects(model))
  split_names <- strsplit(coefficient_names, ":", fixed = TRUE)

  matched <- vapply(
    split_names,
    function(parts) {
      length(parts) == length(required_parts) &&
        setequal(parts, required_parts)
    },
    logical(1)
  )

  result <- coefficient_names[matched]

  if (length(result) != 1) {
    stop(
      "Expected exactly one coefficient containing: ",
      paste(required_parts, collapse = " + "),
      ". Found: ",
      paste(result, collapse = ", ")
    )
  }

  result
}

extract_contrast <- function(
    model,
    weights,
    analysis,
    clinical_outcome,
    component,
    exposure,
    effect,
    participants,
    observations,
    wmh_measurements
) {
  beta <- fixed_effects(model)
  variance_covariance <- as.matrix(vcov(model))

  missing_terms <- setdiff(names(weights), names(beta))
  if (length(missing_terms) > 0) {
    stop(
      "Missing coefficient(s): ",
      paste(missing_terms, collapse = ", ")
    )
  }

  contrast <- setNames(rep(0, length(beta)), names(beta))
  contrast[names(weights)] <- weights

  estimate <- sum(contrast * beta)
  variance <- as.numeric(
    t(contrast) %*% variance_covariance %*% contrast
  )
  standard_error <- sqrt(variance)
  statistic <- estimate / standard_error

  if (inherits(model, "lm") && !inherits(model, "merMod")) {
    degrees_freedom <- df.residual(model)
    critical_value <- qt(0.975, df = degrees_freedom)
    p_value <- 2 * pt(
      abs(statistic),
      df = degrees_freedom,
      lower.tail = FALSE
    )
  } else {
    critical_value <- qnorm(0.975)
    p_value <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
  }

  tibble(
    Analysis = analysis,
    Clinical_outcome = clinical_outcome,
    Component = component,
    Exposure = unname(WMH_LABELS[exposure]),
    Exposure_variable = exposure,
    Effect = effect,
    Participants = participants,
    Observations = observations,
    WMH_measurements = wmh_measurements,
    Beta = estimate,
    Standard_error = standard_error,
    CI_lower = estimate - critical_value * standard_error,
    CI_upper = estimate + critical_value * standard_error,
    P_value = p_value
  )
}

make_sample_row <- function(
    data,
    analysis,
    clinical_outcome,
    component,
    exposure,
    observations_are_wmh_measurements,
    candidate_data = NULL,
    minimum_observations_per_participant = 1L
) {
  if (is.null(candidate_data)) {
    candidate_data <- data
  }

  participant_observation_counts <- data %>%
    count(SubjID, name = "Participant_observations")

  group_counts <- data %>%
    distinct(SubjID, Group) %>%
    count(Group, name = "Participants")

  nonconverter_n <- group_counts %>%
    filter(Group == "Non-converter") %>%
    pull(Participants)

  phenoconverter_n <- group_counts %>%
    filter(Group == "Phenoconverter") %>%
    pull(Participants)

  if (length(nonconverter_n) == 0) nonconverter_n <- 0L
  if (length(phenoconverter_n) == 0) phenoconverter_n <- 0L

  tibble(
    Analysis = analysis,
    Clinical_outcome = clinical_outcome,
    Component = component,
    Exposure = unname(WMH_LABELS[exposure]),
    Exposure_variable = exposure,
    Candidate_participants = n_distinct(candidate_data$SubjID),
    Candidate_observations = nrow(candidate_data),
    Minimum_observations_required =
      minimum_observations_per_participant,
    Participants_excluded_for_minimum =
      n_distinct(candidate_data$SubjID) - n_distinct(data$SubjID),
    Participants = n_distinct(data$SubjID),
    Observations = nrow(data),
    WMH_measurements = if (
      observations_are_wmh_measurements
    ) nrow(data) else n_distinct(data$SubjID),
    Nonconverter_participants = nonconverter_n,
    Phenoconverter_participants = phenoconverter_n,
    Participants_with_1_observation = sum(
      participant_observation_counts$Participant_observations == 1
    ),
    Participants_with_2plus_observations = sum(
      participant_observation_counts$Participant_observations >= 2
    )
  )
}

check_two_groups <- function(data, model_name) {
  observed_groups <- unique(as.character(data$Group))
  observed_groups <- observed_groups[!is.na(observed_groups)]

  if (length(observed_groups) != 2) {
    stop(
      model_name,
      " does not contain both phenoconversion groups after complete-case ",
      "filtering. Observed: ",
      paste(observed_groups, collapse = ", ")
    )
  }
}


# ============================================================================
# 3. VALIDATE THE EXISTING OBJECTS
# ============================================================================

if (!exists("prodromal_long")) {
  stop("prodromal_long does not exist. Run your data-preparation code first.")
}

if (!exists("fpca_results")) {
  stop("fpca_results does not exist. Load your completed FPCA objects first.")
}

prodromal_long <- prodromal_model_data

colnames(prodromal_long)[7] <- "EVENT_ID"


required_source_variables <- c(
  "SubjID",
  "EVENT_ID",
  "Age",
  "years_from_baseline",
  "Sx",
  "Education_Years",
  "Equivalent_Vascular_RF",
  "Native_v",
  unname(WMH_SOURCE_COLUMNS),
  unname(CLINICAL_OUTCOMES)
)

group_source_variable <- if (
  "phenoconversion_group" %in% names(prodromal_long)
) {
  "phenoconversion_group"
} else if (
  "phenoconversion" %in% names(prodromal_long)
) {
  "phenoconversion"
} else {
  stop(
    "prodromal_long must contain phenoconversion_group or phenoconversion."
  )
}

assert_columns(
  prodromal_long,
  c(required_source_variables, group_source_variable),
  "prodromal_long"
)

missing_fpca_domains <- setdiff(
  names(FPCA_SPECIFICATIONS),
  names(fpca_results)
)

if (length(missing_fpca_domains) > 0) {
  stop(
    "fpca_results is missing domains: ",
    paste(missing_fpca_domains, collapse = ", ")
  )
}


# ============================================================================
# 4. BUILD ONE CLEAN VISIT-LEVEL SOURCE TABLE
# ============================================================================

group_input <- prodromal_long[[group_source_variable]]

analysis_long_pre_group <- prodromal_long %>%
  transmute(
    SubjID = as.character(SubjID),
    EVENT_ID = toupper(trimws(as.character(EVENT_ID))),
    Age = numeric_clean(Age),
    supplied_time = numeric_clean(years_from_baseline),
    Group_at_row = normalise_group(group_input),
    Sx = Sx,
    Education_Years = numeric_clean(Education_Years),
    Equivalent_Vascular_RF = Equivalent_Vascular_RF,
    Native_v_visit = numeric_clean(Native_v) / 1000,
    WMH_total = safe_log1p(WMH_native_roi),
    pWMH = safe_log1p(WMH_peri_roi),
    dWMH = safe_log1p(WMH_deep_roi),
    MDS_UPDRS3_OFF_SCORE = numeric_clean(MDS_UPDRS3_OFF_SCORE),
    scopa = numeric_clean(scopa),
    MoCA = numeric_clean(MoCA),
    Cognition_EF = numeric_clean(Cognition_EF),
    Cognition_M = numeric_clean(Cognition_M)
  ) %>%
  filter(
    !is.na(SubjID),
    SubjID != "",
    !is.na(EVENT_ID),
    EVENT_ID != "",
    is.finite(Age),
    Age >= MIN_AGE,
    Age <= MAX_AGE
  )

# Stop on duplicated participant-visit keys before any modelling join.
assert_unique_keys(
  analysis_long_pre_group,
  c("SubjID", "EVENT_ID"),
  "analysis_long_pre_group"
)

# Phenoconversion group is an ever-converted participant-level classification.
# If any row identifies a participant as a phenoconverter, that participant is
# classified as a phenoconverter for all models.
participant_group <- analysis_long_pre_group %>%
  group_by(SubjID) %>%
  summarise(
    Group = case_when(
      any(Group_at_row == "Phenoconverter", na.rm = TRUE) ~
        "Phenoconverter",
      any(Group_at_row == "Non-converter", na.rm = TRUE) ~
        "Non-converter",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Non-converter", "Phenoconverter")
    )
  )

assert_unique_keys(participant_group, "SubjID", "participant_group")

analysis_long <- analysis_long_pre_group %>%
  select(-Group_at_row) %>%
  left_join(participant_group, by = "SubjID")

assert_no_row_expansion(
  nrow(analysis_long_pre_group),
  nrow(analysis_long),
  "participant-group join"
)

assert_unique_keys(
  analysis_long,
  c("SubjID", "EVENT_ID"),
  "analysis_long"
)


# ============================================================================
# 5. STRICT BASELINE DATA: EVENT_ID == BL
# ============================================================================

baseline_clinical_data <- analysis_long %>%
  filter(EVENT_ID == BASELINE_EVENT) %>%
  transmute(
    SubjID,
    Group,
    baseline_age = Age,
    Sx,
    Education_Years,
    Equivalent_Vascular_RF,
    baseline_updrs3 = MDS_UPDRS3_OFF_SCORE,
    baseline_scopa = scopa,
    baseline_moca = MoCA,
    baseline_EF = Cognition_EF,
    baseline_M = Cognition_M
  )

assert_unique_keys(
  baseline_clinical_data,
  "SubjID",
  "baseline_clinical_data"
)

baseline_wmh_data <- analysis_long %>%
  filter(EVENT_ID == BASELINE_EVENT) %>%
  transmute(
    SubjID,
    WMH_total,
    pWMH,
    dWMH,
    Native_v_baseline = Native_v_visit
  )

assert_unique_keys(
  baseline_wmh_data,
  "SubjID",
  "baseline_wmh_data"
)

baseline_subject_data <- baseline_clinical_data %>%
  inner_join(baseline_wmh_data, by = "SubjID")

assert_unique_keys(
  baseline_subject_data,
  "SubjID",
  "baseline_subject_data"
)


# ============================================================================
# 6. DEFINE TIME AND CREATE AUDITED MODEL SOURCES
# ============================================================================

# Baseline age comes only from the verified BL row. Supplied follow-up time is
# used where available; otherwise Age - baseline_age is used.
analysis_long_with_time <- analysis_long %>%
  left_join(
    baseline_clinical_data %>% select(SubjID, baseline_age),
    by = "SubjID"
  ) %>%
  mutate(
    Time = ifelse(
      is.finite(supplied_time),
      supplied_time,
      Age - baseline_age
    )
  ) %>%
  select(-supplied_time)

assert_no_row_expansion(
  nrow(analysis_long),
  nrow(analysis_long_with_time),
  "baseline-age join"
)

assert_unique_keys(
  analysis_long_with_time,
  c("SubjID", "EVENT_ID"),
  "analysis_long_with_time"
)

# Clinical visits contain every available assessment; there is no >=2 filter.
clinical_visit_data <- analysis_long_with_time %>%
  select(
    SubjID,
    EVENT_ID,
    Time,
    all_of(unname(CLINICAL_OUTCOMES))
  )

assert_unique_keys(
  clinical_visit_data,
  c("SubjID", "EVENT_ID"),
  "clinical_visit_data"
)

# Visit WMH data are separate so that the longitudinal join is visibly and
# explicitly performed using SubjID + EVENT_ID.
wmh_visit_data <- analysis_long_with_time %>%
  select(
    SubjID,
    EVENT_ID,
    Native_v_visit,
    all_of(names(WMH_LABELS))
  )

assert_unique_keys(
  wmh_visit_data,
  c("SubjID", "EVENT_ID"),
  "wmh_visit_data"
)

# Panel B source: one unique baseline WMH row is intentionally copied over the
# participant's repeated clinical assessments by a many-to-one SubjID join.
baseline_lmm_source <- clinical_visit_data %>%
  inner_join(baseline_subject_data, by = "SubjID")

assert_no_row_expansion(
  nrow(clinical_visit_data),
  nrow(baseline_lmm_source),
  "baseline-WMH LMM join"
)

assert_unique_keys(
  baseline_lmm_source,
  c("SubjID", "EVENT_ID"),
  "baseline_lmm_source"
)

# Panel C source: clinical and WMH values must have the exact same visit key.
longitudinal_lmm_source <- clinical_visit_data %>%
  inner_join(
    wmh_visit_data,
    by = c("SubjID", "EVENT_ID")
  )

assert_no_row_expansion(
  nrow(clinical_visit_data),
  nrow(longitudinal_lmm_source),
  "visit-matched clinical-WMH join"
)

longitudinal_rows_before_covariate_join <- nrow(longitudinal_lmm_source)

longitudinal_lmm_source <- longitudinal_lmm_source %>%
  inner_join(
    baseline_clinical_data %>%
      select(SubjID, Group, all_of(MODEL_ADJUSTMENTS)),
    by = "SubjID"
  )

assert_no_row_expansion(
  longitudinal_rows_before_covariate_join,
  nrow(longitudinal_lmm_source),
  "longitudinal baseline-covariate join"
)

assert_unique_keys(
  longitudinal_lmm_source,
  c("SubjID", "EVENT_ID"),
  "longitudinal_lmm_source"
)


# ============================================================================
# 7. PANEL A: EXISTING FPCA SCORES AND BASELINE WMH INTERACTIONS
# ============================================================================

fpca_models <- list()
fpca_result_rows <- list()
fpca_sample_rows <- list()
result_index <- 1L
sample_index <- 1L

for (domain in names(FPCA_SPECIFICATIONS)) {
  specification <- FPCA_SPECIFICATIONS[[domain]]
  fpca_models[[domain]] <- list()

  assert_columns(
    fpca_results[[domain]]$scores,
    c("SubjID", specification$fpc1, specification$fpc2),
    paste0("fpca_results$", domain, "$scores")
  )

  for (exposure in names(WMH_LABELS)) {
    fpca_models[[domain]][[exposure]] <- list()

    for (component in c("FPC1", "FPC2")) {
      score_variable <- specification[[tolower(component)]]

      score_data <- fpca_results[[domain]]$scores %>%
        transmute(
          SubjID = as.character(SubjID),
          FPC_score = numeric_clean(.data[[score_variable]])
        )

      assert_unique_keys(
        score_data,
        "SubjID",
        paste(domain, component, "FPCA scores")
      )

      fpca_rows_before_join <- nrow(score_data)

      joined_data <- score_data %>%
        inner_join(baseline_subject_data, by = "SubjID")

      assert_no_row_expansion(
        fpca_rows_before_join,
        nrow(joined_data),
        paste(domain, component, exposure, "FPCA baseline join")
      )

      assert_unique_keys(
        joined_data,
        "SubjID",
        paste(domain, component, exposure, "FPCA model source")
      )

      required_variables <- unique(c(
        "SubjID",
        "Group",
        "FPC_score",
        exposure,
        specification$baseline,
        MODEL_ADJUSTMENTS,
        "Native_v_baseline"
      ))

      model_data <- complete_model_data(
        joined_data,
        required_variables
      ) %>%
        mutate(
          Group = factor(
            Group,
            levels = c("Non-converter", "Phenoconverter")
          )
        )

      check_two_groups(
        model_data,
        paste(domain, component, exposure, "FPCA model")
      )

      model_formula <- reformulate(
        termlabels = unique(c(
          paste0("Group * ", exposure),
          specification$baseline,
          MODEL_ADJUSTMENTS,
          "Native_v_baseline"
        )),
        response = "FPC_score"
      )

      model <- lm(
        model_formula,
        data = model_data,
        na.action = na.fail
      )

      if (nobs(model) != nrow(model_data)) {
        stop("FPCA model sample count mismatch.")
      }

      fpca_models[[domain]][[exposure]][[component]] <- model

      interaction_term <- find_exact_term(
        model,
        c("GroupPhenoconverter", exposure)
      )

      n_participants <- n_distinct(model_data$SubjID)

      fpca_result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, exposure),
        analysis = "FPCA: baseline WMH",
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = component,
        exposure = exposure,
        effect = "WMH association: non-converters",
        participants = n_participants,
        observations = NA_integer_,
        wmh_measurements = n_participants
      )
      result_index <- result_index + 1L

      fpca_result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, interaction_term),
        analysis = "FPCA: baseline WMH",
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = component,
        exposure = exposure,
        effect = "Interaction: phenoconverters vs non-converters",
        participants = n_participants,
        observations = NA_integer_,
        wmh_measurements = n_participants
      )
      result_index <- result_index + 1L

      fpca_result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(
          c(1, 1),
          c(exposure, interaction_term)
        ),
        analysis = "FPCA: baseline WMH",
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = component,
        exposure = exposure,
        effect = "WMH association: phenoconverters",
        participants = n_participants,
        observations = NA_integer_,
        wmh_measurements = n_participants
      )
      result_index <- result_index + 1L

      fpca_sample_rows[[sample_index]] <- make_sample_row(
        data = model_data,
        analysis = "FPCA: baseline WMH",
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = component,
        exposure = exposure,
        observations_are_wmh_measurements = FALSE
      )
      sample_index <- sample_index + 1L
    }
  }
}

fpca_baseline_wmh_results <- bind_rows(fpca_result_rows)
fpca_model_sample_sizes <- bind_rows(fpca_sample_rows)


# ============================================================================
# 8. MIXED-MODEL FITTER USED FOR PANELS B AND C
#
# The fitted fixed-effects structure is:
#   Time * Group * WMH
#
# This estimates:
#   1. WMH association at Time = 0 in non-converters
#   2. Group x WMH difference at Time = 0
#   3. WMH association at Time = 0 in phenoconverters
#   4. Time x WMH association in non-converters
#   5. Time x Group x WMH interaction
#   6. Time x WMH association in phenoconverters
#
# The three-way interaction tests whether the association between WMH burden
# and longitudinal clinical change differs between phenoconverters and
# non-converters.
#
# No within-person/between-person decomposition is performed.
# ============================================================================

fit_lmm_panel <- function(
    model_source,
    analysis_label,
    native_volume_variable,
    observations_are_wmh_measurements,
    minimum_observations_per_participant = 1L
) {
  if (
    length(minimum_observations_per_participant) != 1L ||
      is.na(minimum_observations_per_participant) ||
      minimum_observations_per_participant < 1L
  ) {
    stop("minimum_observations_per_participant must be at least 1.")
  }

  model_store <- list()
  result_rows <- list()
  sample_rows <- list()
  result_index <- 1L
  sample_index <- 1L

  for (domain in names(CLINICAL_OUTCOMES)) {
    outcome <- CLINICAL_OUTCOMES[[domain]]
    model_store[[domain]] <- list()

    for (exposure in names(WMH_LABELS)) {
      required_variables <- unique(c(
        "SubjID",
        "Group",
        "Time",
        outcome,
        exposure,
        MODEL_ADJUSTMENTS,
        native_volume_variable
      ))

      candidate_model_data <- complete_model_data(
        model_source,
        required_variables
      ) %>%
        mutate(
          SubjID = factor(SubjID),
          Group = factor(
            Group,
            levels = c("Non-converter", "Phenoconverter")
          )
        )

      # Apply the visit-count rule only after outcome/exposure-specific
      # complete-case filtering. Therefore, every retained participant has the
      # required number of observations actually used by this exact model.
      model_data <- candidate_model_data %>%
        group_by(SubjID) %>%
        filter(n() >= minimum_observations_per_participant) %>%
        ungroup() %>%
        droplevels()

      if (nrow(model_data) == 0) {
        stop(
          "No eligible complete cases for ",
          domain,
          " / ",
          exposure,
          " after requiring at least ",
          minimum_observations_per_participant,
          " observation(s) per participant."
        )
      }

      check_two_groups(
        model_data,
        paste(analysis_label, domain, exposure)
      )

      model_formula <- reformulate(
        termlabels = unique(c(
          paste0("Time * Group * ", exposure),
          MODEL_ADJUSTMENTS,
          native_volume_variable,
          "(1 | SubjID)"
        )),
        response = outcome
      )

      model <- lmer(
        model_formula,
        data = model_data,
        REML = FALSE,
        na.action = na.fail,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )

      if (nobs(model) != nrow(model_data)) {
        stop(
          analysis_label,
          " sample count mismatch for ",
          domain,
          " / ",
          exposure,
          "."
        )
      }

      model_frame_used <- model.frame(model)
      model_participants <- n_distinct(model_frame_used$SubjID)
      model_observations <- nrow(model_frame_used)
      model_wmh_measurements <- if (
        observations_are_wmh_measurements
      ) model_observations else model_participants

      model_store[[domain]][[exposure]] <- model

      # --------------------------------------------------------------
      # Identify the two-way and three-way interaction coefficients.
      # Because Non-converter is the reference group:
      #
      # exposure
      #   = WMH association at Time = 0 in non-converters
      #
      # GroupPhenoconverter:exposure
      #   = difference in the Time = 0 WMH association between groups
      #
      # Time:exposure
      #   = effect of WMH on annual clinical change in non-converters
      #
      # Time:GroupPhenoconverter:exposure
      #   = difference between groups in the WMH-associated annual change
      # --------------------------------------------------------------

      group_wmh_term <- find_exact_term(
        model,
        c("GroupPhenoconverter", exposure)
      )

      time_wmh_term <- find_exact_term(
        model,
        c("Time", exposure)
      )

      time_group_wmh_term <- find_exact_term(
        model,
        c("Time", "GroupPhenoconverter", exposure)
      )

      # --------------------------------------------------------------
      # A. WMH associations at Time = 0
      # --------------------------------------------------------------

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, exposure),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "WMH association at Time=0: non-converters",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, group_wmh_term),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "Group x WMH interaction at Time=0",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(
          c(1, 1),
          c(exposure, group_wmh_term)
        ),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "WMH association at Time=0: phenoconverters",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      # --------------------------------------------------------------
      # B. WMH associations with longitudinal change
      #
      # These are the key progression terms.
      # A positive beta means that greater WMH is associated with a more
      # positive annual change in the clinical outcome.
      # A negative beta means that greater WMH is associated with a more
      # negative annual change in the clinical outcome.
      #
      # Interpretation still depends on outcome direction:
      #   - UPDRS/SCOPA: positive change = worsening
      #   - MoCA/EF/Memory: negative change = worsening
      # --------------------------------------------------------------

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, time_wmh_term),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "Time x WMH: non-converters",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(1, time_group_wmh_term),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "Time x Group x WMH interaction",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      result_rows[[result_index]] <- extract_contrast(
        model = model,
        weights = setNames(
          c(1, 1),
          c(time_wmh_term, time_group_wmh_term)
        ),
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        effect = "Time x WMH: phenoconverters",
        participants = model_participants,
        observations = model_observations,
        wmh_measurements = model_wmh_measurements
      )
      result_index <- result_index + 1L

      sample_rows[[sample_index]] <- make_sample_row(
        data = model_data,
        analysis = analysis_label,
        clinical_outcome = OUTCOME_LABELS[[domain]],
        component = NA_character_,
        exposure = exposure,
        observations_are_wmh_measurements =
          observations_are_wmh_measurements,
        candidate_data = candidate_model_data,
        minimum_observations_per_participant =
          minimum_observations_per_participant
      )
      sample_index <- sample_index + 1L
    }
  }

  list(
    models = model_store,
    results = bind_rows(result_rows),
    sample_sizes = bind_rows(sample_rows)
  )
}


# ============================================================================
# 9. PANEL B: REPEATED CLINICAL OUTCOMES WITH BASELINE WMH
# ============================================================================

baseline_lmm_output <- fit_lmm_panel(
  model_source = baseline_lmm_source,
  analysis_label = "LMM: baseline WMH",
  native_volume_variable = "Native_v_baseline",
  observations_are_wmh_measurements = FALSE,
  minimum_observations_per_participant = 1L
)

baseline_lmm_models <- baseline_lmm_output$models
baseline_wmh_lmm_results <- baseline_lmm_output$results
baseline_lmm_sample_sizes <- baseline_lmm_output$sample_sizes


# ============================================================================
# 10. PANEL C: REPEATED CLINICAL OUTCOMES WITH VISIT-MATCHED WMH
#
# Eligibility is model-specific: each participant must contribute at least two
# complete, visit-matched WMH-clinical observations to that exact model.
# ============================================================================

longitudinal_lmm_output <- fit_lmm_panel(
  model_source = longitudinal_lmm_source,
  analysis_label = "LMM: longitudinal visit-matched WMH (2+ visits)",
  native_volume_variable = "Native_v_visit",
  observations_are_wmh_measurements = TRUE,
  minimum_observations_per_participant = 2L
)

longitudinal_lmm_models <- longitudinal_lmm_output$models
longitudinal_wmh_lmm_results <- longitudinal_lmm_output$results
longitudinal_lmm_sample_sizes <- longitudinal_lmm_output$sample_sizes


# ============================================================================
# 11. FDR-ADJUSTED P VALUES AND FORMATTED RESULTS TABLES
#
# FDR is calculated separately for each analysis panel and reported effect.
# ============================================================================

add_fdr <- function(data) {
  data %>%
    group_by(Analysis, Effect) %>%
    mutate(
      FDR_adjusted_P = p.adjust(P_value, method = "BH")
    ) %>%
    ungroup()
}

fpca_baseline_wmh_results <- add_fdr(fpca_baseline_wmh_results)
baseline_wmh_lmm_results <- add_fdr(baseline_wmh_lmm_results)
longitudinal_wmh_lmm_results <- add_fdr(longitudinal_wmh_lmm_results)

all_wmh_results_numeric <- bind_rows(
  fpca_baseline_wmh_results,
  baseline_wmh_lmm_results,
  longitudinal_wmh_lmm_results
)

format_results <- function(data) {
  data %>%
    mutate(
      Component = replace_na(Component, ""),
      Sample = case_when(
        is.na(Observations) ~ paste0(
          Participants,
          " participants; ",
          WMH_measurements,
          " baseline WMH measurements"
        ),
        Analysis == "LMM: baseline WMH" ~ paste0(
          Participants,
          " participants; ",
          Observations,
          " clinical observations; ",
          WMH_measurements,
          " baseline WMH measurements"
        ),
        TRUE ~ paste0(
          Participants,
          " participants; ",
          Observations,
          " matched clinical-WMH observations"
        )
      ),
      `beta (95% CI)` = sprintf(
        "%.3f (%.3f to %.3f)",
        Beta,
        CI_lower,
        CI_upper
      ),
      `P value` = format_p(P_value),
      `FDR-adjusted P value` = format_p(FDR_adjusted_P)
    ) %>%
    select(
      Analysis,
      `Clinical outcome` = Clinical_outcome,
      Component,
      Exposure,
      Effect,
      Sample,
      `beta (95% CI)`,
      `P value`,
      `FDR-adjusted P value`
    )
}

fpca_baseline_wmh_table <- format_results(fpca_baseline_wmh_results)
baseline_wmh_lmm_table <- format_results(baseline_wmh_lmm_results)
longitudinal_wmh_lmm_table <- format_results(
  longitudinal_wmh_lmm_results
)
comprehensive_wmh_table <- format_results(all_wmh_results_numeric)

# Progression-specific tables: these contain only terms involving Time x WMH.
# They are useful when the primary question is whether WMH is associated with
# clinical progression and whether that association differs by conversion group.
baseline_wmh_progression_results <- baseline_wmh_lmm_results %>%
  filter(grepl("^Time x ", Effect))

longitudinal_wmh_progression_results <- longitudinal_wmh_lmm_results %>%
  filter(grepl("^Time x ", Effect))

baseline_wmh_progression_table <- format_results(
  baseline_wmh_progression_results
)

longitudinal_wmh_progression_table <- format_results(
  longitudinal_wmh_progression_results
)


# ============================================================================
# 12. ACCURATE MODEL SAMPLE-SIZE TABLE
# ============================================================================

model_sample_sizes <- bind_rows(
  fpca_model_sample_sizes,
  baseline_lmm_sample_sizes,
  longitudinal_lmm_sample_sizes
) %>%
  arrange(
    Analysis,
    Clinical_outcome,
    Component,
    Exposure
  )


# ============================================================================
# 13. BASELINE AND LONGITUDINAL WMH DESCRIPTIVE SUMMARIES
# ============================================================================

summarise_wmh <- function(
    data,
    include_observations,
    minimum_observations_per_participant = 1L
) {
  result <- data %>%
    select(SubjID, Group, all_of(names(WMH_LABELS))) %>%
    pivot_longer(
      cols = all_of(names(WMH_LABELS)),
      names_to = "Exposure_variable",
      values_to = "WMH"
    ) %>%
    filter(!is.na(WMH)) %>%
    group_by(SubjID, Group, Exposure_variable) %>%
    filter(n() >= minimum_observations_per_participant) %>%
    ungroup() %>%
    group_by(Group, Exposure_variable) %>%
    summarise(
      Participants = n_distinct(SubjID),
      Observations = n(),
      Median = median(WMH),
      Q1 = quantile(WMH, 0.25),
      Q3 = quantile(WMH, 0.75),
      .groups = "drop"
    ) %>%
    mutate(
      Exposure = unname(WMH_LABELS[Exposure_variable]),
      `Median (IQR)` = sprintf(
        "%.3f (%.3f to %.3f)",
        Median,
        Q1,
        Q3
      )
    )

  if (!include_observations) {
    result <- result %>% select(-Observations)
  }

  result %>%
    select(
      Group,
      Exposure,
      Participants,
      any_of("Observations"),
      `Median (IQR)`
    )
}

baseline_wmh_descriptive_table <- baseline_subject_data %>%
  summarise_wmh(
    include_observations = FALSE,
    minimum_observations_per_participant = 1L
  )

longitudinal_wmh_descriptive_table <- longitudinal_lmm_source %>%
  summarise_wmh(
    include_observations = TRUE,
    minimum_observations_per_participant = 2L
  )


# ============================================================================
# 14. JOIN AND LONGITUDINAL-WMH AUDIT TABLES
# ============================================================================

join_audit <- tibble(
  Dataset = c(
    "Age-restricted source visits",
    "Baseline clinical rows (BL)",
    "Baseline WMH rows (BL)",
    "Baseline subject rows",
    "All clinical assessment rows",
    "Baseline-WMH LMM source rows",
    "Visit-matched longitudinal LMM source rows"
  ),
  Rows = c(
    nrow(analysis_long),
    nrow(baseline_clinical_data),
    nrow(baseline_wmh_data),
    nrow(baseline_subject_data),
    nrow(clinical_visit_data),
    nrow(baseline_lmm_source),
    nrow(longitudinal_lmm_source)
  ),
  Participants = c(
    n_distinct(analysis_long$SubjID),
    n_distinct(baseline_clinical_data$SubjID),
    n_distinct(baseline_wmh_data$SubjID),
    n_distinct(baseline_subject_data$SubjID),
    n_distinct(clinical_visit_data$SubjID),
    n_distinct(baseline_lmm_source$SubjID),
    n_distinct(longitudinal_lmm_source$SubjID)
  ),
  Key = c(
    "SubjID + EVENT_ID",
    "SubjID",
    "SubjID",
    "SubjID",
    "SubjID + EVENT_ID",
    "SubjID + EVENT_ID",
    "SubjID + EVENT_ID"
  )
)

wmh_longitudinal_qc <- map_dfr(
  names(WMH_LABELS),
  function(exposure) {
    longitudinal_lmm_source %>%
      group_by(SubjID) %>%
      summarise(
        Nonmissing_WMH_visits = sum(!is.na(.data[[exposure]])),
        Distinct_WMH_values = n_distinct(
          .data[[exposure]],
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      summarise(
        Exposure = unname(WMH_LABELS[exposure]),
        Participants_with_any_WMH = sum(Nonmissing_WMH_visits >= 1),
        Participants_with_2plus_WMH_visits = sum(
          Nonmissing_WMH_visits >= 2
        ),
        Participants_with_2plus_distinct_WMH_values = sum(
          Distinct_WMH_values >= 2
        )
      )
  }
)


# ============================================================================
# 15. VIEW RESULTS AND EXAMPLE MODELS
# ============================================================================

print(join_audit)
print(wmh_longitudinal_qc)
print(model_sample_sizes)
print(fpca_baseline_wmh_table)
print(baseline_wmh_lmm_table)
print(longitudinal_wmh_lmm_table)

# Key progression results:
print(baseline_wmh_progression_table)
print(longitudinal_wmh_progression_table)

print(baseline_wmh_descriptive_table)
print(longitudinal_wmh_descriptive_table)

# Example detailed model summaries:
summary(fpca_models$UPDRS$WMH_total$FPC1)

# Examples showing the full Time x Group x WMH fixed-effects structure:
summary(baseline_lmm_models$UPDRS$WMH_total)
summary(baseline_lmm_models$UPDRS$pWMH)
summary(baseline_lmm_models$UPDRS$dWMH)

summary(longitudinal_lmm_models$UPDRS$WMH_total)
summary(longitudinal_lmm_models$UPDRS$pWMH)
summary(longitudinal_lmm_models$UPDRS$dWMH)


# ============================================================================
# 16. SAVE ALL TABLES
# ============================================================================

write.csv(
  fpca_baseline_wmh_table,
  "FPCA_baseline_WMH_interactions_all_assessments.csv",
  row.names = FALSE
)

write.csv(
  baseline_wmh_lmm_table,
  "LMM_baseline_WMH_interactions_all_assessments.csv",
  row.names = FALSE
)

write.csv(
  longitudinal_wmh_lmm_table,
  "LMM_longitudinal_WMH_interactions_all_assessments.csv",
  row.names = FALSE
)

write.csv(
  baseline_wmh_progression_table,
  "LMM_baseline_WMH_TIME_interactions.csv",
  row.names = FALSE
)

write.csv(
  longitudinal_wmh_progression_table,
  "LMM_longitudinal_WMH_TIME_interactions.csv",
  row.names = FALSE
)

write.csv(
  comprehensive_wmh_table,
  "Comprehensive_global_WMH_results_formatted.csv",
  row.names = FALSE
)

write.csv(
  all_wmh_results_numeric,
  "Comprehensive_global_WMH_results_numeric.csv",
  row.names = FALSE
)

write.csv(
  model_sample_sizes,
  "WMH_model_sample_sizes.csv",
  row.names = FALSE
)

write.csv(
  baseline_wmh_descriptive_table,
  "Baseline_WMH_descriptive_summary.csv",
  row.names = FALSE
)

write.csv(
  longitudinal_wmh_descriptive_table,
  "Longitudinal_WMH_descriptive_summary.csv",
  row.names = FALSE
)

write.csv(
  join_audit,
  "WMH_join_audit.csv",
  row.names = FALSE
)

write.csv(
  wmh_longitudinal_qc,
  "WMH_longitudinal_QC_counts.csv",
  row.names = FALSE
)


## one more try:

# ============================================================================
# 17. ANNUALIZED WMH CHANGE AND RELATIONSHIP WITH FPCA SCORES
#
# Paste this section AFTER the previous global WMH / FPCA script.
#
# PRIMARY DEFINITION OF ANNUALIZED WMH CHANGE:
#   For each participant and each WMH compartment:
#
#       raw WMH volume ~ Time
#
#   The participant-specific ordinary least-squares slope is used as the
#   annualized WMH change.
#
#   Therefore:
#     positive slope = increasing WMH burden over time
#     negative slope = decreasing WMH burden over time
#
#   The slope uses ALL available WMH visits, not only first and last scans.
#
# ELIGIBILITY:
#   - At least 2 non-missing WMH measurements
#   - At least 2 distinct time points
#   - Positive follow-up duration
#
# FPCA MODELS:
#
#   FPC score ~ Group * annualized_WMH_change
#               + baseline WMH
#               + baseline clinical outcome
#               + baseline age
#               + sex
#               + education
#               + vascular risk
#               + baseline native brain volume
#
# Total, periventricular and deep WMH are modelled separately.
#
# The annualized change is calculated from the ORIGINAL, NON-LOG WMH volumes.
# Baseline WMH adjustment uses the existing log(1 + WMH) baseline variable,
# consistent with the previous analysis.
# ============================================================================


# ============================================================================
# 17.1 SETTINGS
# ============================================================================

ANNUAL_WMH_SPECS <- list(
  WMH_total = list(
    raw_variable = "WMH_native_roi",
    label = "Total WMH"
  ),
  pWMH = list(
    raw_variable = "WMH_peri_roi",
    label = "Periventricular WMH"
  ),
  dWMH = list(
    raw_variable = "WMH_deep_roi",
    label = "Deep WMH"
  )
)


# ============================================================================
# 17.2 CHECK THAT REQUIRED OBJECTS FROM PREVIOUS SCRIPT EXIST
# ============================================================================

required_existing_objects <- c(
  "prodromal_long",
  "analysis_long_with_time",
  "baseline_subject_data",
  "fpca_results",
  "FPCA_SPECIFICATIONS",
  "OUTCOME_LABELS",
  "MODEL_ADJUSTMENTS"
)

missing_existing_objects <- required_existing_objects[
  !vapply(required_existing_objects, exists, logical(1))
]

if (length(missing_existing_objects) > 0) {
  stop(
    "Run the previous WMH/FPCA script first. Missing object(s): ",
    paste(missing_existing_objects, collapse = ", ")
  )
}

required_existing_functions <- c(
  "numeric_clean",
  "assert_unique_keys",
  "assert_no_row_expansion",
  "check_two_groups",
  "complete_model_data",
  "fixed_effects",
  "find_exact_term",
  "format_p"
)

missing_existing_functions <- required_existing_functions[
  !vapply(required_existing_functions, exists, logical(1))
]

if (length(missing_existing_functions) > 0) {
  stop(
    "Run the previous WMH/FPCA script first. Missing function(s): ",
    paste(missing_existing_functions, collapse = ", ")
  )
}


# ============================================================================
# 17.3 BUILD RAW VISIT-LEVEL WMH TABLE
#
# The earlier analysis_long table contains log(1 + WMH).
# Here we return to the original WMH volume columns so annualized change
# remains directly interpretable as volume change per year.
# ============================================================================

raw_wmh_visit_source <- prodromal_long %>%
  transmute(
    SubjID = as.character(SubjID),
    EVENT_ID = toupper(trimws(as.character(EVENT_ID))),
    Age_check = numeric_clean(Age),

    WMH_total_raw = numeric_clean(WMH_native_roi),
    pWMH_raw = numeric_clean(WMH_peri_roi),
    dWMH_raw = numeric_clean(WMH_deep_roi)
  ) %>%
  filter(
    !is.na(SubjID),
    SubjID != "",
    !is.na(EVENT_ID),
    EVENT_ID != "",
    is.finite(Age_check),
    Age_check >= MIN_AGE,
    Age_check <= MAX_AGE
  )

assert_unique_keys(
  raw_wmh_visit_source,
  c("SubjID", "EVENT_ID"),
  "raw_wmh_visit_source"
)


# Add the exact Time variable already used in the longitudinal models.

rows_before_time_join <- nrow(raw_wmh_visit_source)

raw_wmh_visit_source <- raw_wmh_visit_source %>%
  inner_join(
    analysis_long_with_time %>%
      select(
        SubjID,
        EVENT_ID,
        Time
      ),
    by = c("SubjID", "EVENT_ID")
  )

assert_no_row_expansion(
  rows_before_time_join,
  nrow(raw_wmh_visit_source),
  "raw WMH + Time join"
)

assert_unique_keys(
  raw_wmh_visit_source,
  c("SubjID", "EVENT_ID"),
  "raw_wmh_visit_source with Time"
)


# ============================================================================
# 17.4 FUNCTION TO CALCULATE PARTICIPANT-SPECIFIC ANNUALIZED WMH CHANGE
#
# Annualized change = OLS slope:
#
#       WMH = intercept + beta(Time)
#
# beta is the annualized WMH change.
#
# For participants with exactly 2 visits, this is identical to:
#
#       (WMH2 - WMH1) / (Time2 - Time1)
#
# For participants with 3+ visits, all visits contribute to the slope.
# ============================================================================

calculate_annualized_wmh <- function(
    data,
    value_variable,
    prefix
) {

  result <- data %>%
    transmute(
      SubjID,
      Time = numeric_clean(Time),
      WMH_value = numeric_clean(.data[[value_variable]])
    ) %>%
    filter(
      is.finite(Time),
      is.finite(WMH_value),
      WMH_value >= 0
    ) %>%
    arrange(
      SubjID,
      Time
    ) %>%
    group_by(SubjID) %>%
    summarise(

      n_visits = n(),

      n_distinct_timepoints = n_distinct(Time),

      first_time = min(Time),

      last_time = max(Time),

      followup_years =
        max(Time) - min(Time),

      first_WMH =
        WMH_value[which.min(Time)[1]],

      last_WMH =
        WMH_value[which.max(Time)[1]],

      # OLS slope using every available WMH measurement.
      annualized_change = if (
        n_distinct(Time) >= 2 &&
        is.finite(var(Time)) &&
        var(Time) > 0
      ) {
        cov(Time, WMH_value) / var(Time)
      } else {
        NA_real_
      },

      # Also save first-to-last rate for QC/comparison.
      first_last_annualized_change = if (
        max(Time) > min(Time)
      ) {
        (
          WMH_value[which.max(Time)[1]] -
          WMH_value[which.min(Time)[1]]
        ) /
          (
            max(Time) - min(Time)
          )
      } else {
        NA_real_
      },

      .groups = "drop"
    ) %>%

    # Require true longitudinal WMH information.
    filter(
      n_visits >= 2,
      n_distinct_timepoints >= 2,
      followup_years > 0,
      is.finite(annualized_change)
    )

  # Give every output variable an exposure-specific name.
  names_to_rename <- setdiff(
    names(result),
    "SubjID"
  )

  names(result)[match(
    names_to_rename,
    names(result)
  )] <- paste0(
    prefix,
    "_",
    names_to_rename
  )

  result
}


# ============================================================================
# 17.5 CALCULATE TOTAL, PERIVENTRICULAR AND DEEP WMH ANNUALIZED CHANGE
# ============================================================================

annualized_total_wmh <- calculate_annualized_wmh(
  data = raw_wmh_visit_source,
  value_variable = "WMH_total_raw",
  prefix = "WMH_total"
)

annualized_periventricular_wmh <- calculate_annualized_wmh(
  data = raw_wmh_visit_source,
  value_variable = "pWMH_raw",
  prefix = "pWMH"
)

annualized_deep_wmh <- calculate_annualized_wmh(
  data = raw_wmh_visit_source,
  value_variable = "dWMH_raw",
  prefix = "dWMH"
)


# Combine into one participant-level table.

annualized_wmh_change_subject_level <- annualized_total_wmh %>%
  full_join(
    annualized_periventricular_wmh,
    by = "SubjID"
  ) %>%
  full_join(
    annualized_deep_wmh,
    by = "SubjID"
  ) %>%
  arrange(SubjID)

assert_unique_keys(
  annualized_wmh_change_subject_level,
  "SubjID",
  "annualized_wmh_change_subject_level"
)


# ============================================================================
# 17.6 ADD PHENOCONVERSION GROUP FOR DESCRIPTIVE/QC OUTPUT
# ============================================================================

annualized_wmh_change_subject_level_with_group <-
  annualized_wmh_change_subject_level %>%
  left_join(
    baseline_subject_data %>%
      select(
        SubjID,
        Group
      ),
    by = "SubjID"
  )

assert_unique_keys(
  annualized_wmh_change_subject_level_with_group,
  "SubjID",
  "annualized_wmh_change_subject_level_with_group"
)


# ============================================================================
# 17.7 BASIC QC COUNTS
# ============================================================================

annualized_wmh_qc <- tibble(
  Exposure = c(
    "Total WMH",
    "Periventricular WMH",
    "Deep WMH"
  ),

  Participants_with_2plus_WMH_visits = c(
    sum(
      !is.na(
        annualized_wmh_change_subject_level$
          WMH_total_annualized_change
      )
    ),
    sum(
      !is.na(
        annualized_wmh_change_subject_level$
          pWMH_annualized_change
      )
    ),
    sum(
      !is.na(
        annualized_wmh_change_subject_level$
          dWMH_annualized_change
      )
    )
  ),

  Median_number_of_WMH_visits = c(
    median(
      annualized_wmh_change_subject_level$
        WMH_total_n_visits,
      na.rm = TRUE
    ),
    median(
      annualized_wmh_change_subject_level$
        pWMH_n_visits,
      na.rm = TRUE
    ),
    median(
      annualized_wmh_change_subject_level$
        dWMH_n_visits,
      na.rm = TRUE
    )
  ),

  Median_WMH_followup_years = c(
    median(
      annualized_wmh_change_subject_level$
        WMH_total_followup_years,
      na.rm = TRUE
    ),
    median(
      annualized_wmh_change_subject_level$
        pWMH_followup_years,
      na.rm = TRUE
    ),
    median(
      annualized_wmh_change_subject_level$
        dWMH_followup_years,
      na.rm = TRUE
    )
  )
)


# ============================================================================
# 17.8 DESCRIPTIVE SUMMARY OF ANNUALIZED WMH CHANGE
# ============================================================================

make_annualized_descriptive <- function(
    data,
    exposure,
    label
) {

  slope_variable <- paste0(
    exposure,
    "_annualized_change"
  )

  followup_variable <- paste0(
    exposure,
    "_followup_years"
  )

  visit_variable <- paste0(
    exposure,
    "_n_visits"
  )

  data %>%
    select(
      SubjID,
      Group,
      all_of(slope_variable),
      all_of(followup_variable),
      all_of(visit_variable)
    ) %>%
    rename(
      Annualized_change =
        all_of(slope_variable),

      Followup_years =
        all_of(followup_variable),

      WMH_visits =
        all_of(visit_variable)
    ) %>%
    filter(
      !is.na(Group),
      is.finite(Annualized_change)
    ) %>%
    group_by(Group) %>%
    summarise(
      Participants = n_distinct(SubjID),

      Median_annualized_change =
        median(
          Annualized_change,
          na.rm = TRUE
        ),

      Q1_annualized_change =
        quantile(
          Annualized_change,
          0.25,
          na.rm = TRUE
        ),

      Q3_annualized_change =
        quantile(
          Annualized_change,
          0.75,
          na.rm = TRUE
        ),

      Mean_annualized_change =
        mean(
          Annualized_change,
          na.rm = TRUE
        ),

      SD_annualized_change =
        sd(
          Annualized_change,
          na.rm = TRUE
        ),

      Minimum_annualized_change =
        min(
          Annualized_change,
          na.rm = TRUE
        ),

      Maximum_annualized_change =
        max(
          Annualized_change,
          na.rm = TRUE
        ),

      Median_WMH_visits =
        median(
          WMH_visits,
          na.rm = TRUE
        ),

      Median_followup_years =
        median(
          Followup_years,
          na.rm = TRUE
        ),

      .groups = "drop"
    ) %>%
    mutate(
      Exposure = label,

      `Median annualized change (IQR)` =
        sprintf(
          "%.4f (%.4f to %.4f)",
          Median_annualized_change,
          Q1_annualized_change,
          Q3_annualized_change
        )
    ) %>%
    select(
      Exposure,
      Group,
      Participants,
      `Median annualized change (IQR)`,
      Mean_annualized_change,
      SD_annualized_change,
      Median_WMH_visits,
      Median_followup_years,
      Minimum_annualized_change,
      Maximum_annualized_change
    )
}


annualized_wmh_descriptive_table <- bind_rows(

  make_annualized_descriptive(
    annualized_wmh_change_subject_level_with_group,
    exposure = "WMH_total",
    label = "Total WMH"
  ),

  make_annualized_descriptive(
    annualized_wmh_change_subject_level_with_group,
    exposure = "pWMH",
    label = "Periventricular WMH"
  ),

  make_annualized_descriptive(
    annualized_wmh_change_subject_level_with_group,
    exposure = "dWMH",
    label = "Deep WMH"
  )
)


# ============================================================================
# 17.9 CONTRAST EXTRACTION FOR ANNUALIZED-CHANGE FPCA MODELS
# ============================================================================

extract_annualized_contrast <- function(
    model,
    weights,
    clinical_outcome,
    component,
    exposure,
    exposure_label,
    effect,
    participants
) {

  beta <- fixed_effects(model)

  variance_covariance <-
    as.matrix(vcov(model))

  missing_terms <- setdiff(
    names(weights),
    names(beta)
  )

  if (length(missing_terms) > 0) {
    stop(
      "Missing coefficient(s): ",
      paste(
        missing_terms,
        collapse = ", "
      )
    )
  }

  contrast <- setNames(
    rep(
      0,
      length(beta)
    ),
    names(beta)
  )

  contrast[names(weights)] <-
    weights

  estimate <-
    sum(
      contrast * beta
    )

  variance <-
    as.numeric(
      t(contrast) %*%
        variance_covariance %*%
        contrast
    )

  standard_error <-
    sqrt(variance)

  statistic <-
    estimate /
    standard_error

  degrees_freedom <-
    df.residual(model)

  critical_value <-
    qt(
      0.975,
      df = degrees_freedom
    )

  p_value <-
    2 *
    pt(
      abs(statistic),
      df = degrees_freedom,
      lower.tail = FALSE
    )

  tibble(
    Analysis =
      "FPCA: annualized WMH change",

    Clinical_outcome =
      clinical_outcome,

    Component =
      component,

    Exposure =
      exposure_label,

    Exposure_variable =
      exposure,

    Effect =
      effect,

    Participants =
      participants,

    Beta =
      estimate,

    Standard_error =
      standard_error,

    CI_lower =
      estimate -
      critical_value *
      standard_error,

    CI_upper =
      estimate +
      critical_value *
      standard_error,

    P_value =
      p_value
  )
}


# ============================================================================
# 17.10 FIT FPC ~ ANNUALIZED WMH CHANGE * PHENOCONVERSION MODELS
#
# IMPORTANT:
#
# We adjust for BASELINE WMH burden of the same compartment.
#
# Example:
#
# FPC1_z ~ Group * annualized_total_WMH_change
#          + baseline_total_WMH
#          + baseline_UPDRS
#          + baseline_age
#          + Sx
#          + Education_Years
#          + Equivalent_Vascular_RF
#          + Native_v_baseline
#
# This asks whether WMH ACCUMULATION is associated with FPC score above and
# beyond the participant's starting WMH burden.
# ============================================================================

annualized_fpca_models <- list()

annualized_fpca_result_rows <- list()

annualized_fpca_sample_rows <- list()

annualized_result_index <- 1L

annualized_sample_index <- 1L


for (domain in names(FPCA_SPECIFICATIONS)) {

  specification <-
    FPCA_SPECIFICATIONS[[domain]]

  annualized_fpca_models[[domain]] <-
    list()


  # --------------------------------------------------------------------------
  # Get this clinical domain's FPC scores.
  # --------------------------------------------------------------------------

  score_variables <- c(
    specification$fpc1,
    specification$fpc2
  )

  assert_columns(
    fpca_results[[domain]]$scores,
    c(
      "SubjID",
      score_variables
    ),
    paste0(
      "fpca_results$",
      domain,
      "$scores"
    )
  )


  for (exposure in names(ANNUAL_WMH_SPECS)) {

    annualized_fpca_models[[domain]][[exposure]] <-
      list()

    exposure_label <-
      ANNUAL_WMH_SPECS[[exposure]]$label

    annualized_change_variable <-
      paste0(
        exposure,
        "_annualized_change"
      )

    followup_variable <-
      paste0(
        exposure,
        "_followup_years"
      )

    visit_count_variable <-
      paste0(
        exposure,
        "_n_visits"
      )


    for (component in c(
      "FPC1",
      "FPC2"
    )) {

      score_variable <-
        specification[[
          tolower(component)
        ]]


      # ----------------------------------------------------------------------
      # One row per subject containing the selected FPC score.
      # ----------------------------------------------------------------------

      score_data <-
        fpca_results[[domain]]$scores %>%
        transmute(
          SubjID =
            as.character(SubjID),

          FPC_score =
            numeric_clean(
              .data[[score_variable]]
            )
        )

      assert_unique_keys(
        score_data,
        "SubjID",
        paste(
          domain,
          component,
          "annualized-WMH FPCA scores"
        )
      )


      # ----------------------------------------------------------------------
      # Join:
      #   FPC score
      #   baseline data
      #   annualized WMH slope
      # ----------------------------------------------------------------------

      joined_data <-
        score_data %>%
        inner_join(
          baseline_subject_data,
          by = "SubjID"
        ) %>%
        inner_join(
          annualized_wmh_change_subject_level %>%
            select(
              SubjID,
              all_of(
                annualized_change_variable
              ),
              all_of(
                followup_variable
              ),
              all_of(
                visit_count_variable
              )
            ),
          by = "SubjID"
        )


      assert_unique_keys(
        joined_data,
        "SubjID",
        paste(
          domain,
          component,
          exposure,
          "annualized WMH model source"
        )
      )


      # ----------------------------------------------------------------------
      # Variables required for this exact model.
      # ----------------------------------------------------------------------

      required_variables <-
        unique(
          c(
            "SubjID",
            "Group",
            "FPC_score",

            # Main exposure:
            annualized_change_variable,

            # Adjust for starting WMH burden:
            exposure,

            # Baseline score of matching clinical outcome:
            specification$baseline,

            MODEL_ADJUSTMENTS,

            "Native_v_baseline",

            followup_variable,

            visit_count_variable
          )
        )


      model_data <-
        complete_model_data(
          joined_data,
          required_variables
        ) %>%
        mutate(
          Group = factor(
            Group,
            levels = c(
              "Non-converter",
              "Phenoconverter"
            )
          )
        )


      check_two_groups(
        model_data,
        paste(
          domain,
          component,
          exposure,
          "annualized WMH FPCA model"
        )
      )


      # ----------------------------------------------------------------------
      # Model.
      # ----------------------------------------------------------------------

      model_formula <-
        reformulate(
          termlabels =
            unique(
              c(
                paste0(
                  "Group * ",
                  annualized_change_variable
                ),

                # Baseline burden of same WMH compartment.
                exposure,

                specification$baseline,

                MODEL_ADJUSTMENTS,

                "Native_v_baseline"
              )
            ),

          response =
            "FPC_score"
        )


      model <-
        lm(
          model_formula,
          data = model_data,
          na.action = na.fail
        )


      if (
        nobs(model) !=
        nrow(model_data)
      ) {
        stop(
          "Annualized WMH FPCA model sample-count mismatch for ",
          domain,
          " / ",
          component,
          " / ",
          exposure
        )
      }


      annualized_fpca_models[[
        domain
      ]][[
        exposure
      ]][[
        component
      ]] <- model


      # ----------------------------------------------------------------------
      # Identify Group x annualized WMH interaction.
      # ----------------------------------------------------------------------

      interaction_term <-
        find_exact_term(
          model,
          c(
            "GroupPhenoconverter",
            annualized_change_variable
          )
        )


      n_participants <-
        n_distinct(
          model_data$SubjID
        )


      # ----------------------------------------------------------------------
      # 1. Association in NON-CONVERTERS
      #
      # Reference group = Non-converters.
      # Therefore coefficient of annualized WMH change alone is the effect in
      # non-converters.
      # ----------------------------------------------------------------------

      annualized_fpca_result_rows[[
        annualized_result_index
      ]] <-
        extract_annualized_contrast(
          model = model,

          weights =
            setNames(
              1,
              annualized_change_variable
            ),

          clinical_outcome =
            OUTCOME_LABELS[[domain]],

          component =
            component,

          exposure =
            annualized_change_variable,

          exposure_label =
            exposure_label,

          effect =
            "Annualized WMH change association: non-converters",

          participants =
            n_participants
        )

      annualized_result_index <-
        annualized_result_index + 1L


      # ----------------------------------------------------------------------
      # 2. PHENOCONVERSION INTERACTION
      #
      # Does the annualized WMH-change association differ between
      # phenoconverters and non-converters?
      # ----------------------------------------------------------------------

      annualized_fpca_result_rows[[
        annualized_result_index
      ]] <-
        extract_annualized_contrast(
          model = model,

          weights =
            setNames(
              1,
              interaction_term
            ),

          clinical_outcome =
            OUTCOME_LABELS[[domain]],

          component =
            component,

          exposure =
            annualized_change_variable,

          exposure_label =
            exposure_label,

          effect =
            "Interaction: phenoconverters vs non-converters",

          participants =
            n_participants
        )

      annualized_result_index <-
        annualized_result_index + 1L


      # ----------------------------------------------------------------------
      # 3. Association in PHENOCONVERTERS
      #
      # Phenoconverter effect:
      #
      # annualized-change main effect + interaction
      # ----------------------------------------------------------------------

      annualized_fpca_result_rows[[
        annualized_result_index
      ]] <-
        extract_annualized_contrast(
          model = model,

          weights =
            setNames(
              c(
                1,
                1
              ),
              c(
                annualized_change_variable,
                interaction_term
              )
            ),

          clinical_outcome =
            OUTCOME_LABELS[[domain]],

          component =
            component,

          exposure =
            annualized_change_variable,

          exposure_label =
            exposure_label,

          effect =
            "Annualized WMH change association: phenoconverters",

          participants =
            n_participants
        )

      annualized_result_index <-
        annualized_result_index + 1L


      # ----------------------------------------------------------------------
      # Save accurate sample information.
      # ----------------------------------------------------------------------

      group_counts <-
        model_data %>%
        distinct(
          SubjID,
          Group
        ) %>%
        count(
          Group,
          name = "Participants"
        )


      nonconverter_n <-
        group_counts %>%
        filter(
          Group ==
            "Non-converter"
        ) %>%
        pull(
          Participants
        )

      phenoconverter_n <-
        group_counts %>%
        filter(
          Group ==
            "Phenoconverter"
        ) %>%
        pull(
          Participants
        )


      if (
        length(nonconverter_n) == 0
      ) {
        nonconverter_n <- 0L
      }

      if (
        length(phenoconverter_n) == 0
      ) {
        phenoconverter_n <- 0L
      }


      annualized_fpca_sample_rows[[
        annualized_sample_index
      ]] <-
        tibble(
          Analysis =
            "FPCA: annualized WMH change",

          Clinical_outcome =
            OUTCOME_LABELS[[domain]],

          Component =
            component,

          Exposure =
            exposure_label,

          Participants =
            n_participants,

          Nonconverter_participants =
            nonconverter_n,

          Phenoconverter_participants =
            phenoconverter_n,

          Median_WMH_visits =
            median(
              model_data[[
                visit_count_variable
              ]],
              na.rm = TRUE
            ),

          Minimum_WMH_visits =
            min(
              model_data[[
                visit_count_variable
              ]],
              na.rm = TRUE
            ),

          Maximum_WMH_visits =
            max(
              model_data[[
                visit_count_variable
              ]],
              na.rm = TRUE
            ),

          Median_WMH_followup_years =
            median(
              model_data[[
                followup_variable
              ]],
              na.rm = TRUE
            ),

          Minimum_WMH_followup_years =
            min(
              model_data[[
                followup_variable
              ]],
              na.rm = TRUE
            ),

          Maximum_WMH_followup_years =
            max(
              model_data[[
                followup_variable
              ]],
              na.rm = TRUE
            )
        )

      annualized_sample_index <-
        annualized_sample_index + 1L
    }
  }
}


# ============================================================================
# 17.11 COMBINE MODEL RESULTS
# ============================================================================

fpca_annualized_wmh_results <-
  bind_rows(
    annualized_fpca_result_rows
  )

fpca_annualized_wmh_sample_sizes <-
  bind_rows(
    annualized_fpca_sample_rows
  )


# ============================================================================
# 17.12 FDR CORRECTION
#
# Same strategy as your previous script:
# FDR is calculated separately for each reported Effect.
#
# This means the following are corrected separately:
#
#   1. Annualized change association in non-converters
#   2. Phenoconverter vs non-converter interaction
#   3. Annualized change association in phenoconverters
# ============================================================================

fpca_annualized_wmh_results <-
  fpca_annualized_wmh_results %>%
  group_by(
    Analysis,
    Effect
  ) %>%
  mutate(
    FDR_adjusted_P =
      p.adjust(
        P_value,
        method = "BH"
      )
  ) %>%
  ungroup()


# ============================================================================
# 17.13 FORMATTED RESULTS TABLE
# ============================================================================

fpca_annualized_wmh_table <-
  fpca_annualized_wmh_results %>%
  mutate(

    `beta (95% CI)` =
      sprintf(
        "%.3f (%.3f to %.3f)",
        Beta,
        CI_lower,
        CI_upper
      ),

    `P value` =
      format_p(
        P_value
      ),

    `FDR-adjusted P value` =
      format_p(
        FDR_adjusted_P
      ),

    Sample =
      paste0(
        Participants,
        " participants"
      )
  ) %>%
  select(
    Analysis,

    `Clinical outcome` =
      Clinical_outcome,

    Component,

    Exposure,

    Effect,

    Sample,

    `beta (95% CI)`,

    `P value`,

    `FDR-adjusted P value`
  )


# ============================================================================
# 17.14 OPTIONAL: TABLE SHOWING ONLY THE PHENOCONVERSION INTERACTION
#
# This is particularly useful for your forest plots / dissertation summary.
# ============================================================================

fpca_annualized_wmh_interaction_table <-
  fpca_annualized_wmh_results %>%
  filter(
    Effect ==
      "Interaction: phenoconverters vs non-converters"
  ) %>%
  mutate(

    `beta (95% CI)` =
      sprintf(
        "%.3f (%.3f to %.3f)",
        Beta,
        CI_lower,
        CI_upper
      ),

    `P value` =
      format_p(
        P_value
      ),

    `FDR-adjusted P value` =
      format_p(
        FDR_adjusted_P
      )
  ) %>%
  select(
    `Clinical outcome` =
      Clinical_outcome,

    Component,

    Exposure,

    Participants,

    `beta (95% CI)`,

    `P value`,

    `FDR-adjusted P value`
  )


# ============================================================================
# 17.15 PRINT RESULTS
# ============================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "ANNUALIZED WMH CHANGE QC\n"
)

cat(
  "============================================================\n\n"
)

print(
  annualized_wmh_qc
)


cat(
  "\n\n============================================================\n"
)

cat(
  "ANNUALIZED WMH CHANGE DESCRIPTIVES\n"
)

cat(
  "============================================================\n\n"
)

print(
  annualized_wmh_descriptive_table
)


cat(
  "\n\n============================================================\n"
)

cat(
  "FPCA ~ ANNUALIZED WMH CHANGE RESULTS\n"
)

cat(
  "============================================================\n\n"
)

print(
  fpca_annualized_wmh_table
)


cat(
  "\n\n============================================================\n"
)

cat(
  "PHENOCONVERSION INTERACTIONS ONLY\n"
)

cat(
  "============================================================\n\n"
)

print(
  fpca_annualized_wmh_interaction_table
)


cat(
  "\n\n============================================================\n"
)

cat(
  "MODEL SAMPLE SIZES\n"
)

cat(
  "============================================================\n\n"
)

print(
  fpca_annualized_wmh_sample_sizes
)


# ============================================================================
# 17.16 EXAMPLE DETAILED MODEL SUMMARIES
# ============================================================================

# UPDRS FPC1 ~ annualized TOTAL WMH change
summary(
  annualized_fpca_models$
    UPDRS$
    WMH_total$
    FPC1
)

# UPDRS FPC2 ~ annualized PERIVENTRICULAR WMH change
summary(
  annualized_fpca_models$
    UPDRS$
    pWMH$
    FPC2
)

# Executive-function FPC1 ~ annualized DEEP WMH change
summary(
  annualized_fpca_models$
    EF$
    dWMH$
    FPC1
)


# ============================================================================
# 17.17 SAVE EVERYTHING
# ============================================================================

# Participant-level annualized change data.
write.csv(
  annualized_wmh_change_subject_level_with_group,
  "Annualized_WMH_change_subject_level.csv",
  row.names = FALSE
)


# Basic QC counts.
write.csv(
  annualized_wmh_qc,
  "Annualized_WMH_change_QC.csv",
  row.names = FALSE
)


# Descriptive statistics by phenoconversion group.
write.csv(
  annualized_wmh_descriptive_table,
  "Annualized_WMH_change_descriptive_summary.csv",
  row.names = FALSE
)


# Full numeric model results.
write.csv(
  fpca_annualized_wmh_results,
  "FPCA_annualized_WMH_change_models_numeric.csv",
  row.names = FALSE
)


# Formatted manuscript-style results.
write.csv(
  fpca_annualized_wmh_table,
  "FPCA_annualized_WMH_change_models_formatted.csv",
  row.names = FALSE
)


# Interaction-only table.
write.csv(
  fpca_annualized_wmh_interaction_table,
  "FPCA_annualized_WMH_change_interactions_only.csv",
  row.names = FALSE
)


# Accurate model sample sizes.
write.csv(
  fpca_annualized_wmh_sample_sizes,
  "FPCA_annualized_WMH_change_model_sample_sizes.csv",
  row.names = FALSE
)


cat(
  "\nAnnualized WMH analysis completed successfully.\n"
)

cat(
  "Files written:\n",
  "  Annualized_WMH_change_subject_level.csv\n",
  "  Annualized_WMH_change_QC.csv\n",
  "  Annualized_WMH_change_descriptive_summary.csv\n",
  "  FPCA_annualized_WMH_change_models_numeric.csv\n",
  "  FPCA_annualized_WMH_change_models_formatted.csv\n",
  "  FPCA_annualized_WMH_change_interactions_only.csv\n",
  "  FPCA_annualized_WMH_change_model_sample_sizes.csv\n"
)
