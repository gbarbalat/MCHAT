#this script makes a table of Included vs. Total sample participants (i.e. compares baseline characteristics)
rm(list=ls())
library(dlnm)
library(dplyr)
library(gtsummary)


#Organisation of Tble 1 as in the paper ----
final_orga <- c(
  "EDI",
  "mean_NDVI",
  "M00X4_TAU2010",
  "M00M2_LIEUNAISM",
  "A01M_HxLEARN",
  "A01P_HxLEARN",
  "A01M_HxLANG",
  "A01P_HxLANG",
  "educ_2m",
  "M00M2_CSP1M",
  "M00M2_ENFGANT",
  "M00M2_BMIMAVTG",
  "M00M2_AGEM",
  "M00M2_AGEP",
  "M00M2_FQALCOOL",
  "TOBACCO",
  "M00M3_FQCAFE",
  "M00M3_POISGEN",
  "M00M3_VITB9",
  "M00M3_MGOMEGA3",
  "SEXE_ENF",
  "M02M_TYPALI",
  "A02M_EFVIT",
  "NBLANGMEN_2m",
  "revenu_part_qui_2y",
  "A02M_AGE2A"
)


load("Not_Included_MCHAT.RData")

load("Not_Imputed_MCHAT.RData");Not_Imputed <- Not_Imputed_Outcome
load("Imputed_MCHAT.RData");Imputed <- Imputed_Outcome

exposure_prefix=c("NO2_preN","NO2_postN", "PM2.5_preN", "PM2.5_postN","PM10_preN", "PM10_postN", 
                  "Tmax_preN", "Tmax_postN","Tmin_preN", "Tmin_postN", "Tmean_preN", "Tmean_postN"
)
#transform all numeric into factors, except for exposure and outcomes, also remove exposure data and subject id
Not_Included <- Not_Included %>% select(-all_of(starts_with(exposure_prefix)),id_Dem806_585_GB) %>%
  dplyr::mutate(across(
    .cols = where(is.numeric) & !dplyr::any_of("ALL"),
    .fns = as.factor
  ))
Imputed <- Imputed %>% select(-all_of(starts_with(exposure_prefix)),id_Dem806_585_GB) %>%
  dplyr::mutate(across(
    .cols = where(is.numeric) & !dplyr::any_of("ALL"),
    .fns = as.factor
  )) 
Not_Imputed <- Not_Imputed %>% select(-all_of(starts_with(exposure_prefix)),id_Dem806_585_GB) %>%
  dplyr::mutate(across(
    .cols = where(is.numeric) & !dplyr::any_of("ALL"),
    .fns = as.factor
  ))

Included <- Not_Imputed

#all covariate (col in database)
all_covar <- c(
  
  "EDI",
  "mean_NDVI",
  "NBLANGMEN_2m",
  "M00X4_TAU2010",
  
  #"M00M1_VAGUE",
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  "M00M2_LIEUNAISM",# born in France or another country
  
  "A02M_EFVIT",#efvit = personne avec qui vit l’enfant
  #"M00M2_COUPLE", #live with someone ###
  
  #"M00M2_NIVET", #8 categ
  "educ_2m", ### #higher of meduc_2m feduc_2m,
  
  "M00M2_CSP1M", #9 categ M00M2_PROFESS prefer
  #"M00M2_CSP1P",
  
  "M00M2_ENFGANT", #"M00M2_NBGANT",# Y No DK
  
  "M00M2_BMIMAVTG",
  
  "M00M2_FQALCOOL", #freq etOH (ordinal)
  
  "M00M3_FQCAFE", #6 categ
  
  "M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, always
  
  "M00M3_VITB9", # No Yes DNK
  
  "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Allways
  
  #### 2M
  "M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  "TOBACCO",
  # "M02_EXPTAB", #   "M02M_EXPTAB", "M02P_EXPTAB", #= exposure to tobacco 5 categ
  # "M00M2_TABAG",# N, Y
  # "M00M2_EXPTAB", ### #"M00M2_EXPTABD", "M00M2_EXPTABLF", #passive tobacco (home)--
  
  
  #### 1Y
  #"A01M_HxNDD",
  #"A01P_HxNDD",
  "A01M_HxLANG",
  "A01M_HxLEARN",
  "A01P_HxLANG",
  "A01P_HxLEARN",
  
  #### 2Y
  "A02M_AGE2A",# age in months
  "revenu_part_qui_2y"
  
)

# bind Included with Not_Included to make a Total_Sample
# Add a group indicator to each dataset
Total_Sample <- bind_rows(Not_Imputed, Not_Included)
Included$Group <- "Included"
Total_Sample$Group <- "Total Sample"
Total_Sample <- Total_Sample %>% select(-all_of(starts_with(exposure_prefix)),id_Dem806_585_GB) %>%
  dplyr::mutate(across(
    .cols = where(is.numeric) & !dplyr::any_of("ALL"),
    .fns = as.factor
  ))

# Combine the datasets
combined_data <- bind_rows(Included , Total_Sample)

# Ensure Group is a factor
combined_data$Group <- factor(combined_data$Group, levels = c("Included", "Total Sample"))

# Select only the variables in final_orga and the group indicator
analysis_data <- combined_data %>% select(all_of(final_orga), Group)

# AGE: Replace "29" with "27_28"
analysis_data$A02M_AGE2A[analysis_data$A02M_AGE2A == "29"] <- "27_28"
analysis_data$A02M_AGE2A <- factor(analysis_data$A02M_AGE2A)
# Drop the unused level "29"
analysis_data$A02M_AGE2A <- droplevels(analysis_data$A02M_AGE2A)


# Create the summary table comparing by Group, and add p-values
summary_table <- tbl_summary(
  data = analysis_data,
  by = Group
) %>%
  add_p() # Adds p-values for group comparisons[1][2][3]

# Print in RStudio viewer
summary_table

# To export as HTML
as_gt(summary_table) %>% gt::gtsave("summary_table.html")

#calculate Cramer V to compare included and non-included groups
cramers_v_manual <- function(x, y) {
  tbl <- table(x, y)
  chi2 <- suppressWarnings(chisq.test(tbl, correct = FALSE)$statistic)
  n <- sum(tbl)
  k <- min(nrow(tbl), ncol(tbl))
  v <- sqrt(chi2 / (n * (k - 1)))
  return(as.numeric(v))
}

# Example use:
cramers_v_manual(analysis_data$Group, analysis_data$revenu_part_qui_2y)


factor_vars <- colnames(analysis_data)
cramer_results <- sapply(factor_vars, function(var) {
  cramers_v_manual(analysis_data$Group, analysis_data[[var]])
})

# Format as a named vector or data frame
Cramer <- data.frame(Variable = factor_vars,
                     Cramers_V = round(cramer_results, 3),
                     row.names = NULL)

#extract the gtsummary table into a tibble
tbl_tibble <- summary_table %>%
  as_tibble()

# Join Cramer's V table — match only variable rows (not level rows)
tbl_tibble <- tbl_tibble %>%
  left_join(Cramer, by = c( `**Characteristic**`  = "Variable")) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))


sjPlot::tab_df((tbl_tibble[,"Cramers_V"]))
