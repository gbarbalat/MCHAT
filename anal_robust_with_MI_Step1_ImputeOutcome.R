rm(list=ls())
library(dplyr)
library(dlnm)
library(mice)
library(tidyr)



#mice parameter
m <- 10
maxit <- 10

#Function to better understand the distribution of missing values among variables and participants using naniar
p_miss_indiv <- 30 #if more than p_miss_indiv% of missing value, remove the individual

#function to evaluate missing data
try_missing <- function(data) {
  (miss_var_summary_obs <-naniar::miss_var_summary(data)); print(miss_var_summary_obs)
  (miss_var_table_obs <- naniar::miss_var_table(data) %>%
      mutate(pct_miss_in_var=n_miss_in_var*100/nrow(data), .after=n_miss_in_var))
  print(miss_var_table_obs)
  (miss_case_summary_obs <- naniar::miss_case_summary(data))
  print(miss_case_summary_obs)
  (miss_case_table_obs <- naniar::miss_case_table(data) %>%
      mutate(pct_miss_in_case=n_miss_in_case*100/ncol(data), .after=n_miss_in_case) %>%
      arrange(desc(pct_miss_in_case))
  )
  print(miss_case_table_obs)
  return(miss_case_summary_obs)
}

# Function to calculate LAMBDA and FMI (fraction of missing information)
calculate_fmi <- function(imp, formula) {
  m <- imp$m  # number of imputations
  n <- nrow(imp$data)#sample size
  k <-model.matrix(formula,imp$data) %>% colnames %>% length#nb of param to fit
  
  # Fit model to each imputed dataset
  fits <- lapply(1:m, function(i) {
    data <- complete(imp, i)
    fit <- lm(formula, data = data)
    list(coef = fit$coefficients, se = sqrt(diag(vcov(fit))))
  })
  
  # Extract coefficients and standard errors
  coefs <- sapply(fits, function(x) x$coef)
  ses <- sapply(fits, function(x) x$se)
  
  # Calculate components
  theta_bar <- rowMeans(coefs)
  # Calculate summed differences for each coefficient
  B <- sapply(rownames(coefs), function(coef_name) {
    sum((coefs[coef_name,] - theta_bar[coef_name])^2)/(m-1)
  })
  names(B) <- rownames(coefs)  
  W <- rowMeans(ses^2)#averaged within imputation variance
  
  # Calculate FMI
  T <- W + B + B/m # T <- W+B
  lambda <- (B + B/m)/ T; lambda
  riv <- lambda/(1-lambda); riv
  df_old <- (m-1)/lambda^2; df_old
  df_obs <- (((n-k)+1)/((n-k)+3))*(n-k)*(1-lambda); df_obs
  df_adj <- (df_old*df_obs)/(df_old+df_obs);df_adj
  #The difference between lambda and FMI is that FMI is adjusted for the fact that the number of imputed datasets 
  #that are generated is not unlimitedly large. These measures differ for a small value of the df.
  num_fmi <- riv + 2/(df_adj+3)
  den_fmi <- 1+riv
  fmi <- num_fmi/den_fmi; fmi
  
  return(list(lambda,fmi))
}


#covariate vector
#caution: wave and dad's employment status have been removed (VIF>5)
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
all_num_final <- c("ALL")#numeric value
exposure_prefix <- c("NO2_preN","NO2_postN", "PM2.5_preN", "PM2.5_postN","PM10_preN", "PM10_postN", 
                     "Tmax_preN", "Tmax_postN","Tmin_preN", "Tmin_postN", "Tmean_preN", "Tmean_postN"
)
all_exp_mean <- c("TMM_preN" ,"TMaM_preN","TMiM_preN","TMM_postN", 
                  "TMaM_postN","TMiM_postN","PM2.5M_preN","PM10M_preN","NO2M_preN",         
                  "PM2.5M_postN","PM10M_postN","NO2M_postN")

#load database with missing data
load("df_covar_merged_final_recode.RData")
#load id to be thrown away (multiple births, no QS at 2 years old)
load("some_id.RData")
str(df_covar_merged_final_recode)
#make character as factors
df_covar_merged_final_recode <- df_covar_merged_final_recode %>%
  mutate(across(where(is.character), as.factor))
#take only those with a 2y questionnaire
df_covar_merged_final_recode <- df_covar_merged_final_recode %>%
  filter(!id_Dem806_585_GB %in% id_MB) %>%
  filter(!id_Dem806_585_GB %in% id_2y_NoQS)
str(df_covar_merged_final_recode)
#merge with exposure database
load("MCHAT_exposure.RData")
# 
# load("Not_Imputed_Outcome.RData"); colSums(is.na(Not_Imputed_Outcome)); Outcome_exposure <- Not_Imputed_Outcome
# #load("Not_Imputed_MCHAT.RData")
# # load("Not_Included_MCHAT.RData")
# # 
# # load("df_covar_merged.RData")#18329
# # load("df_covar_merged_easy_recode.RData")#18329
# # load("df_covar_merged_final_recode.RData")#18329
# load("df_covar_merged_final_recode_NA.RData"); str(df_covar_merged_final_recode_NA); summary(df_covar_merged_final_recode_NA$ALL)
# colSums(is.na(df_covar_merged_final_recode_NA))
# load("df_covar_merged_final_recode_NA_Imputed.RData");str(df_covar_merged_final_recode_NA_Imputed); summary(df_covar_merged_final_recode_NA_Imputed$ALL)
# colSums(is.na(df_covar_merged_final_recode_NA_Imputed))


MCHAT_exposure_NO2 <- Outcome_exposure %>%
  select(-garde_alter) %>%
  mutate(NO2_allNA = if_all(starts_with("NO2"), is.na),
         NO2_someNA = {
           NO2_cols <- select(., starts_with("NO2"))
           rowSums(is.na(NO2_cols)) > 0 & rowSums(!is.na(NO2_cols)) > 0
         }
  )

#remove missing exposure data as when missing, it is just too much (whole pre and postnatal periods OR the whole postnatal period)
MCHAT_exposure_NO2 <- MCHAT_exposure_NO2 %>%
  filter(!(NO2_allNA | NO2_someNA))

# merge exposure data with covariate data
merged_ <- df_covar_merged_final_recode %>% select(-ALL) %>%
  inner_join( MCHAT_exposure_NO2 ,by="id_Dem806_585_GB") 


#merged_col_obs select only columns of interest -> alrady done
#merged_gp gp levels so that you have enough N per level -> alrady done
#merged_ignore ignore rows/observations because too many missing values AND perform IPW of staying -> NOT DONE
merged_ignore_withExp <- merged_ %>%
  mutate(TMM_preN=rowMeans(select(merged_,all_of(starts_with("Tmean_preN")))%>% na.omit)) %>%
  mutate(TMaM_preN=rowMeans(select(merged_,all_of(starts_with("Tmax_preN")))%>% na.omit)) %>%
  mutate(TMiM_preN=rowMeans(select(merged_,all_of(starts_with("Tmin_preN")))%>% na.omit)) %>%
  
  mutate(TMM_postN=rowMeans(select(merged_,all_of(starts_with("Tmean_postN")))%>% na.omit)) %>%
  mutate(TMaM_postN=rowMeans(select(merged_,all_of(starts_with("Tmax_postN")))%>% na.omit)) %>%
  mutate(TMiM_postN=rowMeans(select(merged_,all_of(starts_with("Tmin_postN")))%>% na.omit)) %>%
  
  mutate(PM2.5M_preN=rowMeans(select(merged_,all_of(starts_with("PM2.5_preN")))%>% na.omit)) %>%
  mutate(PM10M_preN=rowMeans(select(merged_,all_of(starts_with("PM10_preN")))%>% na.omit)) %>%
  mutate(NO2M_preN=rowMeans(select(merged_,all_of(starts_with("NO2_preN")))%>% na.omit)) %>%
  
  mutate(PM2.5M_postN=rowMeans(select(merged_,all_of(starts_with("PM2.5_postN")))%>% na.omit)) %>%
  mutate(PM10M_postN=rowMeans(select(merged_,all_of(starts_with("PM10_postN")))%>% na.omit)) %>%
  mutate(NO2M_postN=rowMeans(select(merged_,all_of(starts_with("NO2_postN")))%>% na.omit))
  

merged_ignore <- merged_ignore_withExp %>%
  select(-all_of(starts_with(exposure_prefix)), -id_Dem806_585_GB) 

#evaluate missing data
miss_case_summary_obs <- try_missing(merged_ignore)
#cases with less than xxx% missing values
keep_id <- miss_case_summary_obs %>%
  filter(pct_miss<p_miss_indiv) %>%
  dplyr::select(case)
#only select cases with more than xxx% missing values
merged_ignore <- merged_ignore[sort(keep_id$case),] 
merged_ignore_withExp <- merged_ignore_withExp[sort(keep_id$case),] 
fx <- flux(merged_ignore)
plot(fx$influx, fx$outflux, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Influx", ylab = "Outflux", main = "Flux Plot")
text(fx$influx, fx$outflux, row.names(fx), pos = 4, cex = 0.8)

# Replace "29" with "27_28"
merged_ignore$A02M_AGE2A[merged_ignore$A02M_AGE2A == "29"] <- "27_28"

# Drop the unused level "29"
merged_ignore$A02M_AGE2A <- droplevels(merged_ignore$A02M_AGE2A)

#remove cols allNA and someNA
merged_ignore <- merged_ignore %>%
  select(-all_of(ends_with("allNA")),-all_of(ends_with("someNA")))

#imputation step
colSums(is.na(merged_ignore))

set.seed(12345)
merged_imputed <- mice::mice(merged_ignore, m=m, maxit = maxit)
merged_imputed$loggedEvents

formula_FMI <- as.formula(paste0("as.integer(ALL) ~ ",
                                 paste0(all_exp_mean, collapse = "+"),
                                 "+",
                                 paste0(all_covar, collapse = "+"))
)
(FMI_LAMBDA <- calculate_fmi(merged_imputed, formula_FMI))

save.image(file=paste0("Imputed_MCHAT.RData" ))

#end of script