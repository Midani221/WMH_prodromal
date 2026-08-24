library(xlsx)
library(tidyverse)
library(stringr)

#Paths of file containing all datasets
path_data <- "PATH TO SHEETS" # Change to path containing all sheets names
setwd(path_data)

#
variables <- read.xlsx("Variables.xlsx",1)
curated_set <- read.csv("PPMI_Curated_datacut.csv")

PDAQ27_dataset <- read.csv("PDAQ27.csv")
updrs_datasets <- c("MDS-UPDRS_Part_I_PQ","MDS-UPDRS_Part_II_PQ","MDS-UPDRS_Part_III","MDS-UPDRS_Part_IV")

medical_conditions <- read.csv("vascular_risk_factors.csv")
medications <- read.csv("All_Medication_Log.csv")


prodromal_imaging_samples <- read.csv("Prodrome_Dataset_12_22.csv")
imaging_details <- read.csv("Imaging_characteristics.csv")
pd_imaging_samples <- read.csv("PD_dataset_12_22.csv")

###
name_vars <- variables$Variable_PPMI

colns <- which(colnames(curated_set) %in% name_vars)
data <- curated_set[,colns]

# Join updrs to curated set data
updrs_list <- map(
  updrs_datasets,
  ~ read.csv(paste0(.x, ".csv"))
)

# left join all of them to data
data <- reduce(
  updrs_list,
  ~ left_join(.x, .y, by = c("PATNO", "EVENT_ID")),
  .init = data
)

pdaq_score <- PDAQ27_dataset %>%
  select(PATNO, EVENT_ID, PDAQ27)

data <- data %>%
  left_join(pdaq_score, by = c("PATNO", "EVENT_ID"))

prodromal_s <- prodromal_imaging_samples %>%
		select(PATNO,EVENT_ID,T1_FLAIR,T1_FLAIR_DTI,Visit_date_imaging)

pd_s <- pd_imaging_samples %>%
		select(PATNO,EVENT_ID,T1_FLAIR,T1_FLAIR_DTI,Visit_date_imaging)


imaging_s <- add_row(prodromal_s,pd_s)


imaging_details$Field_Strength <- as.numeric(sub(".*Field Strength=([^;]+).*", "\\1", imaging_details$Imaging.Protocol))
imaging_details$Scanner <- sub(".*Manufacturer=([^;]+).*", "\\1", imaging_details$Imaging.Protocol)

imaging_d <- imaging_details %>%
	select(PATNO,Field_Strength,Scanner)

imaging_d <- imaging_d %>%
     distinct(PATNO,.keep_all =TRUE)

imaging_s <- imaging_s %>%
	left_join(imaging_d,by="PATNO")

imaging_s <- imaging_s %>%
     distinct(PATNO,EVENT_ID,.keep_all =TRUE)


data <- data %>%
  left_join(imaging_s, by = c("PATNO", "EVENT_ID"))

data <- data %>%
  left_join(medical_conditions, by = "PATNO")

idx <- match(colnames(data), variables$Variable_PPMI)
new_names <- colnames(data)
new_names[!is.na(idx)] <- variables[[2]][idx[!is.na(idx)]]
colnames(data) <- new_names

df <- as.data.frame(matrix(NA,ncol=length(variables$Variable_New)))
colnames(df) <- variables$Variable_New

df2 <- bind_rows(df,data)
df2 <- df2[,1:length(variables$Variable_New)]
df2 <- df2[-1,]

# Things to fix: cognitive things.

df2$Medicated <- ((df2$LEDD > 0) + df2$COMTi + df2$DA + df2$Amantadine + df2$Anticholinergic + df2$MAOBI) > 0

df2$Cognition_AWM <- ifelse(!is.na(df2$DVT_SDM) & !is.na(df2$DVZ_TMTA),((df2$DVT_SDM-50)/10 + df2$DVZ_TMTA)/2, 
				ifelse(!is.na(df2$DVT_SDM), (df2$DVT_SDM - 50)/10,df2$DVZ_TMTA)) 

df2$Cognition_EF <- ifelse(!is.na(df2$DVT_FAS) & !is.na(df2$DVZ_TMTB),((df2$DVT_FAS-50)/10 + df2$DVZ_TMTB)/2, 
				ifelse(!is.na(df2$DVT_FAS), (df2$DVT_FAS - 50)/10,df2$DVZ_TMTB)) 

df2$Cognition_L <- ifelse(!is.na(df2$DVT_SFTANIM) & !is.na(df2$DVS_BNT),((df2$DVT_SFTANIM-50)/10 + (df2$DVS_BNT-10)/3)/2, 
				ifelse(!is.na(df2$DVT_SFTANIM), (df2$DVT_SFTANIM - 50)/10,(df2$DVS_BNT-10)/3)) 

df2$Cognition_M <- ((df2$DVT_TOTAL_RECALL-50)/10 + (df2$DVT_DELAYED_RECALL-50)/10)/2


df2$Cognition_VF <- (df2$DVS_JLO_MSSAE - 10)/3

# Any vascular RF
vascular_rf_vars <- c(
  "BMI_VascularRF",
  "Smoking",
  "Pre_Borderline_DiabetesMellitus_VascularRF",
  "DiabetesMellitus_VascularRF",
  "Hypertension_VascularRF",
  "Hypercholestrolaemia_VascularRF",
  "Dyslipidaemia_VascularRF",
  "IschaemicHeartDisease_VascularRF",
  "PreviousTIA_Stroke_VascularRF",
  "AtrialFibrillation_VascularRF",
  "HeartDisease_VascularRF",
  "Systemic_Inflammatory_Condition_VascularRF",
  "Cancer_Present"
)

df2 <- df2 %>%
  rowwise() %>%
  mutate(
    Any_VascularRF = case_when(
      any(c_across(all_of(vascular_rf_vars)) == 1, na.rm = TRUE) ~ 1,
      all(is.na(c_across(all_of(vascular_rf_vars)))) ~ NA_real_,
      TRUE ~ 0
    )
  ) %>%
  ungroup()




#Phenoconversion

df2$Visit_date_clinical <- mdy(df2$Visit_date_clinical)

df2 <- df2 %>%
  mutate(
    visit_order = case_when(
      Visit_No == "BL" ~ 0,
      str_detect(Visit_No, "^V\\d+") ~ as.numeric(str_remove(Visit_No, "V")),
      TRUE ~ NA_real_
    )
  )

phenoconversion_summary <- df2 %>%
  filter(Dx== 4) %>%
  arrange(SubjID, visit_order, Visit_date_clinical) %>%
  group_by(SubjID) %>%
  summarise(
    baseline_date = Visit_date_clinical[Visit_No == "BL"][1],
    baseline_PRIMDIAG = PRIM_DIAG[Visit_No == "BL"][1],

    phenoconversion = any(PRIM_DIAG== 1 & Visit_No != "BL", na.rm = TRUE),

    phenoconversion_EVENT_ID = ifelse(
      phenoconversion,
      Visit_No[which(PRIM_DIAG== 1 & Visit_No != "BL")[1]],
      NA_character_
    ),

    phenoconversion_date = ifelse(
      phenoconversion,
      as.character(Visit_date_clinical[which(PRIM_DIAG== 1 & Visit_No != "BL")[1]]),
      NA_character_
    ),

    time_to_phenoconversion_days = ifelse(
      phenoconversion,
      as.numeric(as.Date(phenoconversion_date) - baseline_date),
      NA_real_
    ),

    .groups = "drop"
  ) %>%
  mutate(
    phenoconversion = as.integer(phenoconversion),
    phenoconversion_date = as.Date(phenoconversion_date)
  )

df2 <- df2 %>%
  left_join(
    phenoconversion_summary %>%
      select(
        SubjID,
        phenoconversion,
        phenoconversion_EVENT_ID,
        phenoconversion_date,
        time_to_phenoconversion_days
      ),
    by = "SubjID"
  )


#
write.csv(df2,"Clinical_dataset.csv")


