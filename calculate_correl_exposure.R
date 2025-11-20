#This simply calculates correlations between weekly exposure indicators

rm(list=ls())
library(dplyr)
library(tidyr)

load(file=paste0("Imputed_MCHAT.RData"))

#prenatal and postnatal lags, see define_cx_b_AIC function 
nlags_preN=29; i_lag_max_preN = 30; step_preN=7
nlags_postN=23; i_lag_max_postN = 24; step_postN=30
nlags_postN=19; i_lag_max_postN = 20; step_postN=30 #Discard last 4 months of exposure
nlags_postN=90; i_lag_max_postN = 91; step_postN=7 #weekly 1 year and 3/4


exposure_prefix=c("NO2_preN","NO2_postN", "PM2.5_preN", "PM2.5_postN","PM10_preN", "PM10_postN", 
                  "Tmax_preN", "Tmax_postN","Tmin_preN", "Tmin_postN", "Tmean_preN", "Tmean_postN"
)



# function to Calculate weekly exposure matrix ----
var="Tmin_preN_"; #Tmean_preN_ Tmax_preN_ Tmin_preN_
varfun="ns";nlags=nlags_preN; i_lag_max = i_lag_max_preN; step=step_preN;
weekly_exp_matrix <- function(var,nlags,i_lag_max,step) {
  i_lag=1
  i_d=1
  xx=list()
  while (i_lag<=i_lag_max) {
    x=list()
    for (i_step in 1:step) {
      var_i_d=ifelse(grepl("preN", var, fixed = TRUE),i_d-1+365,i_d-1)
      x[[i_step]]=Imputed_Outcome %>% 
        select(paste0(var,var_i_d))
      i_d=i_d+1
    }
    #xx=do.call(cbind,x)
    xx[[i_lag]]=do.call(cbind,x) %>% rowMeans()
    
    i_lag=i_lag+1
  }
  xxx=do.call(cbind,xx)
}


weekly.NO2_preN <- weekly_exp_matrix(var="NO2_preN_", # NO2_preN_ PM2.5_preN_ PM10_preN_
                               nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.NO2_postN <- weekly_exp_matrix(var="NO2_postN_", # NO2_postN PM2.5_postN PM10_postN
                                nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.NO2 <- cbind(weekly.NO2_preN,weekly.NO2_postN)
  
weekly.PM2.5_preN <- weekly_exp_matrix(var="PM2.5_preN_", 
                                     nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.PM2.5_postN <- weekly_exp_matrix(var="PM2.5_postN_",
                                      nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.PM2.5 <- cbind(weekly.PM2.5_preN,weekly.PM2.5_postN)

weekly.PM10_preN <- weekly_exp_matrix(var="PM10_preN_", 
                                       nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.PM10_postN <- weekly_exp_matrix(var="PM10_postN_",
                                        nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.PM10 <- cbind(weekly.PM10_preN,weekly.PM10_postN)

weekly.Tmax_preN <- weekly_exp_matrix(var="Tmax_preN_", #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.Tmax_postN <- weekly_exp_matrix(var="Tmax_postN_",
                                 nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.Tmax <- cbind(weekly.Tmax_preN,weekly.Tmax_postN)

weekly.Tmin_preN <- weekly_exp_matrix(var="Tmin_preN_", 
                                      nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.Tmin_postN <- weekly_exp_matrix(var="Tmin_postN_",
                                       nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.Tmin <- cbind(weekly.Tmin_preN,weekly.Tmin_postN)

weekly.Tmean_preN <- weekly_exp_matrix(var="Tmean_preN_", 
                                      nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN
)
weekly.Tmean_postN <- weekly_exp_matrix(var="Tmean_postN_",
                                       nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN
)
weekly.Tmean <- cbind(weekly.Tmean_preN,weekly.Tmean_postN)



# concatenate and calculate correlation matrix ----
# List all your weekly pollutant matrices here
exposures <- list(weekly.Tmean, weekly.Tmax, weekly.Tmin,weekly.PM2.5, weekly.PM10, weekly.NO2 )  # Add all your matrices
exposure_names <- c("Tmean", "Tmax", "Tmin", "PM2.5","PM10" , "NO2")  # Corresponding names

# Vectorize each matrix and combine into a data frame
exposure_data <- data.frame(sapply(exposures, function(x) as.vector(x)))
colnames(exposure_data) <- exposure_names

cor_matrix <- cor(exposure_data, use = "complete.obs")

library(sjPlot)
tab_corr(cor_matrix, show.p = FALSE, triangle = "both")

