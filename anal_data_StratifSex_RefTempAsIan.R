rm(list=ls())
library(dlnm)
library(dplyr)
library(tidyr)

#anal with longer prenatal period or not
path <- "" #"" OR "./37WA/"


# STRATIF="SEX" ----
categ="Female"
#name of the pdf file (graphs) to be created
file=paste0("anal_data_Stratif",categ,".pdf")

#df for cx basis
grid1_Poll=NULL#dose-resp relationship for pollution; df=NULL (linear function)
grid2_Poll_preN=3
grid2_Poll_postN=4
grid1_Temp=2#dose-resp relationship for temperature (ns function with INTERCEPT = FALSE)
grid2_Temp=3 #lag-response relationship (ns function with INTERCEPT = TRUE)

nlags_postN=23; i_lag_max_postN = 24; step_postN=30
nlags_postN=19; i_lag_max_postN = 20; step_postN=30 #Discard last 4 months of exposure
nlags_postN=90; i_lag_max_postN = 91; step_postN=7 #weekly 1 year and 3/4

# Ref values (mean and pctiles) ----
#This is the ref values used in cross preds: median and various percentiles
#one for each period (prenat vs. postnat) and each temp variable (Tmax, Tmin, Tmean)
find_ref_AsIan =function(var,i_lag_max,step) {
  i_lag=1
  i_d=1
  xx=list()
  while (i_lag<=i_lag_max) {
    x=list()
    #select 7(preN)/30(postN) days
    for (i_step in 1:step) {
      #caution: conception starts 365 days after initial records 
      var_i_d=ifelse(grepl("preN", var, fixed = TRUE),i_d-1+365,i_d-1)
      x[[i_step]]=Imputed_Outcome %>% 
        select(paste0(var,var_i_d))
      i_d=i_d+1
    }
    #calculate mean of 7(preN)/30(postN) days
    xx[[i_lag]]=do.call(cbind,x) %>% rowMeans()
    
    i_lag=i_lag+1
  }
  #make as long dataframe and calculate quantiles
  xxx=do.call(cbind,xx) %>%
    as.data.frame() %>% pivot_longer(cols=everything(),names_to = "names",values_to = "val") %>%
    reframe(out=quantile(val, probs=c(0,0.01,0.05,0.10,0.20,0.5,0.80, 0.90,0.95,0.99,1)))
  
}

#set ref temperatures and  prenatal DLNM parameters depending on whether anal on longer prenatal periods
#also load appropriate database
if (path=="") {
  load(file=paste0(path,"Imputed_MCHAT.RData"))
  
  #prenatal  lags 
  nlags_preN=29; i_lag_max_preN = 30; step_preN=7
  
  Tmean_preN_pct=find_ref_AsIan("Tmean_preN_",i_lag_max=i_lag_max_preN,step=step_preN)
  Tmin_preN_pct=find_ref_AsIan("Tmin_preN_",i_lag_max=i_lag_max_preN,step=step_preN)
  Tmax_preN_pct=find_ref_AsIan("Tmax_preN_",i_lag_max=i_lag_max_preN,step=step_preN)
  
  Tmean_postN_pct=find_ref_AsIan("Tmean_postN_",i_lag_max=i_lag_max_postN,step=step_postN)
  Tmin_postN_pct=find_ref_AsIan("Tmin_postN_",i_lag_max=i_lag_max_postN,step=step_postN)
  Tmax_postN_pct=find_ref_AsIan("Tmax_postN_",i_lag_max=i_lag_max_postN,step=step_postN)
  
  save(Tmean_preN_pct,Tmin_preN_pct,Tmax_preN_pct,
       Tmean_postN_pct,Tmin_postN_pct,Tmax_postN_pct,
       file="ref_temperature.RData")
  
} else if (path=="./37WA/") {
  
  load( file="ref_temperature.RData")
  load(file=paste0(path,"Imputed_MCHAT_37WA.RData"))
  #inner merging to non-preterm children
  load("id_37SA.RData")
  Imputed_Outcome <- Imputed_Outcome %>%
    filter(id_Dem806_585_GB %in% id_37WA)
  #prenatal  lags 
  nlags_preN=34; i_lag_max_preN = 35; step_preN=7
}

#Filter out by Sex
Imputed_Outcome = Imputed_Outcome %>%
  filter(SEXE_ENF==categ)#Male Female
#filter(mean_NDVI==)#(0.0367,0.393] (0.393,0.509] (0.509,0.897]
#filter(EDI==)#(-1.6,2.06]  (2.06,31.9]  (-9.18,-1.6]


#covariate vector
#caution: wave and dad's employment status have been removed (VIF>5)
all_covar=c(
  
  "EDI",
  "mean_NDVI",
  "NBLANGMEN_2m",
  "M00X4_TAU2010",
  
  #"M00M1_VAGUE",
  
  #"SEXE_ENF",
  
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



exposure_prefix=c("NO2_preN","NO2_postN", "PM2.5_preN", "PM2.5_postN","PM10_preN", "PM10_postN", 
                  "Tmax_preN", "Tmax_postN","Tmin_preN", "Tmin_postN", "Tmean_preN", "Tmean_postN"
)



# Define Cx-b with AIC ----
#first part of function identical to find_ref_AsIan
#second part simply uses crossbasis function from dlnm package
# var="Tmin_preN_"; #Tmean_preN_ Tmax_preN_ Tmin_preN_
# varfun="ns";nlags=nlags_preN; i_lag_max = i_lag_max_preN; step=step_preN;
# grid1=grid1_Temp;grid2=grid2_Temp
define_cx_b_AIC =function(var,varfun,nlags,i_lag_max,step, 
                          grid1, grid2) {
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
  head(xxx)
  colSums(is.na(xxx))
  
  if(varfun=="lin") {argvar=list(fun = varfun)
  } else {argvar=list(fun = varfun, df = grid1)}
  
  
  ###
  sub=xxx[1,];exphist(sub, times=length(sub), lag=length(sub)-1)
  xxx=t(apply(xxx,1, function(sub) exphist(sub, 
                                           times=length(sub),
                                           lag=length(sub)-1)
  )
  )
  xxx[1,]
  ###
  cx_b=crossbasis(x = xxx , 
                  lag = nlags, 
                  argvar = argvar, #ER rel
                  arglag = list(fun = "ns", df = grid2) #LR rel
                  
  )
  #colSums(is.na(basis.Tmax))
}


basis.NO2_preN=define_cx_b_AIC(var="NO2_preN_", # NO2_preN_ PM2.5_preN_ PM10_preN_
                               varfun="lin",
                               nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                               grid1=grid1_Poll,grid2=grid2_Poll_preN
                               
)

basis.NO2_postN=define_cx_b_AIC(var="NO2_postN_", # NO2_postN PM2.5_postN PM10_postN
                                varfun="lin",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                grid1=grid1_Poll,grid2=grid2_Poll_postN
                                
)

basis.PM2.5_preN=define_cx_b_AIC(var="PM2.5_preN_", # NO2_preN_ PM2.5_preN_ PM10_preN_
                                 varfun="lin",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                 grid1=grid1_Poll,grid2=grid2_Poll_preN
                                 
)

basis.PM2.5_postN=define_cx_b_AIC(var="PM2.5_postN_", # NO2_postN PM2.5_postN PM10_postN
                                  varfun="lin",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                  grid1=grid1_Poll,grid2=grid2_Poll_postN
                                  
)

basis.PM10_preN=define_cx_b_AIC(var="PM10_preN_", # NO2_preN_ PM2.5_preN_ PM10_preN_
                                varfun="lin",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                grid1=grid1_Poll,grid2=grid2_Poll_preN
                                
)

basis.PM10_postN=define_cx_b_AIC(var="PM10_postN_", # NO2_postN PM2.5_postN PM10_postN
                                 varfun="lin",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                 grid1=grid1_Poll,grid2=grid2_Poll_postN
                                 
)

basis.Tmax_preN=define_cx_b_AIC(var="Tmax_preN_", #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                varfun="ns",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                grid1=grid1_Temp,grid2=grid2_Temp
                                
)

basis.Tmax_postN=define_cx_b_AIC(var="Tmax_postN_",
                                 varfun="ns",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                 grid1=grid1_Temp,grid2=grid2_Temp
                                 
)

basis.Tmin_preN=define_cx_b_AIC(var="Tmin_preN_", #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                varfun="ns",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                grid1=grid1_Temp,grid2=grid2_Temp
                                
)

basis.Tmin_postN=define_cx_b_AIC(var="Tmin_postN_",
                                 varfun="ns",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                 grid1=grid1_Temp,grid2=grid2_Temp
                                 
)

basis.Tmean_preN=define_cx_b_AIC(var="Tmean_preN_", #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                 varfun="ns",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                 grid1=grid1_Temp,grid2=grid2_Temp
                                 
)

basis.Tmean_postN=define_cx_b_AIC(var="Tmean_postN_",
                                  varfun="ns",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                  grid1=grid1_Temp,grid2=grid2_Temp
                                  
)


#- Make ordered linear ----
# When inspecting the data, we found that some ordered variables might be best considered linear 
Imputed_Outcome$M00M3_POISGEN=factor(Imputed_Outcome$M00M3_POISGEN,
                                   levels=c("0","1","2","3","More_2_W"),
                                   ordered=TRUE)
contrasts(Imputed_Outcome$M00M3_POISGEN,1) <- contr.poly(5)[, 1]

#
Imputed_Outcome$revenu_part_qui_2y=factor(Imputed_Outcome$revenu_part_qui_2y,
                                        levels=c("0","1","2","3","4"),
                                        ordered=TRUE)
contrasts(Imputed_Outcome$revenu_part_qui_2y,1) <- contr.poly(5)[, 1]

#
Imputed_Outcome$educ_2m=factor(Imputed_Outcome$educ_2m,
                             levels=c("_2","3","4","5","6"),
                             ordered=TRUE)
contrasts(Imputed_Outcome$educ_2m,1) <- contr.poly(5)[, 1]


# Formulas ----

### Formula for model - Pollution alone
formula_PM2.5=as.formula(paste0("as.integer(ALL) ~  basis.PM2.5_preN + basis.PM2.5_postN +", 
                                paste0(all_covar, collapse = "+")))

formula_PM10=as.formula(paste0("as.integer(ALL) ~  basis.PM10_preN + basis.PM10_postN +", 
                               paste0(all_covar, collapse = "+")))

formula_NO2=as.formula(paste0("as.integer(ALL) ~  basis.NO2_preN + basis.NO2_postN +", 
                              paste0(all_covar, collapse = "+")))

#Formula for model - Temperature alone
formula_Tmax=as.formula(paste0("as.integer(ALL) ~ 
                           basis.Tmax_preN + basis.Tmax_postN +",
                               paste0(all_covar, collapse = "+")))

formula_Tmin=as.formula(paste0("as.integer(ALL) ~ 
                           basis.Tmin_preN + basis.Tmin_postN +",
                               paste0(all_covar, collapse = "+")))

formula_Tmean=as.formula(paste0("as.integer(ALL) ~ 
                           basis.Tmean_preN + basis.Tmean_postN +",
                                paste0(all_covar, collapse = "+")))

#Formula for model - Pollution and temperature
formula_Tmax_NO2=as.formula(paste0("as.integer(ALL) ~  basis.Tmax_preN + basis.Tmax_postN +", 
                                   "basis.NO2_preN + basis.NO2_postN +",
                                   paste0(all_covar, collapse = "+")))
formula_Tmax_PM2.5=as.formula(paste0("as.integer(ALL) ~  basis.Tmax_preN + basis.Tmax_postN +", 
                                     "basis.PM2.5_preN + basis.PM2.5_postN +",
                                     paste0(all_covar, collapse = "+")))
formula_Tmax_PM10=as.formula(paste0("as.integer(ALL) ~  basis.Tmax_preN + basis.Tmax_postN +", 
                                    "basis.PM10_preN + basis.PM10_postN +",
                                    paste0(all_covar, collapse = "+")))

formula_Tmin_NO2=as.formula(paste0("as.integer(ALL) ~  basis.Tmin_preN + basis.Tmin_postN +", 
                                   "basis.NO2_preN + basis.NO2_postN +",
                                   paste0(all_covar, collapse = "+")))
formula_Tmin_PM2.5=as.formula(paste0("as.integer(ALL) ~  basis.Tmin_preN + basis.Tmin_postN +", 
                                     "basis.PM2.5_preN + basis.PM2.5_postN +",
                                     paste0(all_covar, collapse = "+")))
formula_Tmin_PM10=as.formula(paste0("as.integer(ALL) ~  basis.Tmin_preN + basis.Tmin_postN +", 
                                    "basis.PM10_preN + basis.PM10_postN +",
                                    paste0(all_covar, collapse = "+")))

formula_Tmean_NO2=as.formula(paste0("as.integer(ALL) ~  basis.Tmean_preN + basis.Tmean_postN +", 
                                    "basis.NO2_preN + basis.NO2_postN +",
                                    paste0(all_covar, collapse = "+")))
formula_Tmean_PM2.5=as.formula(paste0("as.integer(ALL) ~  basis.Tmean_preN + basis.Tmean_postN +", 
                                      "basis.PM2.5_preN + basis.PM2.5_postN +",
                                      paste0(all_covar, collapse = "+")))
formula_Tmean_PM10=as.formula(paste0("as.integer(ALL) ~  basis.Tmean_preN + basis.Tmean_postN +", 
                                     "basis.PM10_preN + basis.PM10_postN +",
                                     paste0(all_covar, collapse = "+")))

formula=list(formula_PM2.5,formula_PM10,formula_NO2,
             formula_Tmax,formula_Tmin,formula_Tmean,
             formula_Tmax_NO2,formula_Tmax_PM2.5,formula_Tmax_PM10,
             formula_Tmin_NO2,formula_Tmin_PM2.5,formula_Tmin_PM10,
             formula_Tmean_NO2,formula_Tmean_PM2.5,formula_Tmean_PM10)

# Fit Models ----
mod=list()
for (i in 1:length(formula)) {
  mod[[i]] <- #MASS::glm.nb(#  
    glm(family=quasipoisson(),#poisson gaussian quasipoisson
        formula[[i]],
        na.action = na.omit,
        
        Imputed_Outcome 
    )
  
  summary(mod[[i]])
  print(usefr::AICmodel(mod[[i]]))
  print(car::vif(mod[[i]]))
  
}


# Cross pred ----

### PM2.5 mod 1,8,11,14 (see formula)
cp_PM2.5_preN=function(mod_nb) {crosspred(basis.PM2.5_preN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.PM2.5_preN=lapply(c(1,8,11,14),cp_PM2.5_preN)

cp_PM2.5_postN=function(mod_nb) {crosspred(basis.PM2.5_postN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.PM2.5_postN=lapply(c(1,8,11,14),cp_PM2.5_postN)

### PM10 mod 2,9,12,15 (see formula)
cp_PM10_preN=function(mod_nb) {crosspred(basis.PM10_preN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.PM10_preN=lapply(c(2,9,12,15),cp_PM10_preN)

cp_PM10_postN=function(mod_nb) {crosspred(basis.PM10_postN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.PM10_postN=lapply(c(2,9,12,15),cp_PM10_postN)

### NO2 mod 3,7,10,13 (see formula)
cp_NO2_preN=function(mod_nb) {crosspred(basis.NO2_preN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.NO2_preN=lapply(c(3,7,10,13),cp_NO2_preN)

cp_NO2_postN=function(mod_nb) {crosspred(basis.NO2_postN,mod[[mod_nb]],cen=0,bylag = 1,cumul = TRUE)}
pred.NO2_postN=lapply(c(3,7,10,13),cp_NO2_postN)

# Tmax 4,7,8,9 (see formula)
cp_Tmax_preN=function(mod_nb,cen, at) {crosspred(basis.Tmax_preN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmax_preN=lapply(c(4,7,8,9),cp_Tmax_preN,cen=Tmax_preN_pct[6,], at=Tmax_preN_pct %>% unlist)

cp_Tmax_postN=function(mod_nb,cen, at) {crosspred(basis.Tmax_postN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmax_postN=lapply(c(4,7,8,9),cp_Tmax_postN,cen=Tmax_postN_pct[6,], at=Tmax_postN_pct %>% unlist)

# Tmin 5,10,11,12 (see formula)
cp_Tmin_preN=function(mod_nb,cen,at) {crosspred(basis.Tmin_preN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmin_preN=lapply(c(5,10,11,12),cp_Tmin_preN,cen=Tmin_preN_pct[6,], at=Tmin_preN_pct %>% unlist)

cp_Tmin_postN=function(mod_nb,cen,at) {crosspred(basis.Tmin_postN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmin_postN=lapply(c(5,10,11,12),cp_Tmin_postN,cen=Tmin_postN_pct[6,], at=Tmin_postN_pct %>% unlist)

# Tmean 6,13,14,15 (see formula)
cp_Tmean_preN=function(mod_nb,cen,at) {crosspred(basis.Tmean_preN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmean_preN=lapply(c(6,13,14,15),cp_Tmean_preN,cen=Tmean_preN_pct[6,]%>% unlist,
                       at=Tmean_preN_pct %>% unlist)

cp_Tmean_postN=function(mod_nb,cen,at) {crosspred(basis.Tmean_postN,mod[[mod_nb]],cen=cen,bylag = 1,cumul = TRUE, at=at)}
pred.Tmean_postN=lapply(c(6,13,14,15),cp_Tmean_postN,cen=Tmean_postN_pct[6,] %>% unlist, 
                        at=Tmean_postN_pct %>% unlist)

save.image(file=paste0(categ,".RData"))

# #3D graph
# plot(pred_temp, 
#      xlab = "mean of temp", 
#      ylab = "Lag", 
#      main = "3D graph")
# 
# #Contour plot
# plot(pred.Tmean_preN[[2]],
#      "contour",
#      plot.title = title("Contour plot",
#                         xlab = "mean of temp",
#                         ylab = "Lag"))
# 
# #Cumulative effects across lags
# plot(pred_temp, 
#      "overall",
#      xlab = "mean of temp", 
#      xaxs = "i", 
#      ci.arg = list(col = grey(0.7)), 
#      main = "Cumulative association during preN")

## Time-response curves for specific values ----

plot_pred_time=function(which_cp,cp,name,var,pct=NULL) {
  xlab=ifelse(grepl("preN",name,fixed = TRUE),"Weeks","Months")
  xlim=NULL;
  xlim[1]=if_else(grepl("preN",name,fixed = TRUE),i_lag_max_preN,i_lag_max_postN);
  xlim[2]=0
  main=ifelse(grepl("Tm",name,fixed = TRUE),paste0(name," ",pct,"th pct"),name)
  plot(cp[[which_cp]], 
       "slices", 
       var = var,
       ylab = 'Relative Risk', 
       xlab = xlab, 
       xlim = xlim,#c(38,0)
       xaxs = "i", 
       xaxt = "n", 
       col = "#619CFF", 
       lwd = 2, 
       ci.arg = list(col = "#619CFF40"),
       main = main, 
       col.main = "#619CFF")
}



{
  pdf(file=paste0("Poll_",file))
  
  lapply(1:4,plot_pred_time,pred.PM2.5_preN,"pred.PM2.5_preN",10)
  lapply(1:4,plot_pred_time,pred.PM2.5_postN,"pred.PM2.5_postN",10)
  
  lapply(1:4,plot_pred_time,pred.PM10_preN,"pred.PM10_preN",10)
  lapply(1:4,plot_pred_time,pred.PM10_postN,"pred.PM10_postN",10)
  
  lapply(1:4,plot_pred_time,pred.NO2_preN,"pred.NO2_preN",10)
  lapply(1:4,plot_pred_time,pred.NO2_postN,"pred.NO2_postN",10)
  
  dev.off()
  
  #Tmin
  pdf(file=paste0("Tmin_",file))
  
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmin_preN,"pred.Tmin_preN",unlist(Tmin_preN_pct[10,]), pct="99")
  # 
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmin_postN,"pred.Tmin_postN",unlist(Tmin_postN_pct[10,]), pct="99")
  dev.off()
  
  #Tmax
  pdf(file=paste0("Tmax_",file))
  
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmax_preN,"pred.Tmax_preN",unlist(Tmax_preN_pct[10,]), pct="99")
  # 
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmax_postN,"pred.Tmax_postN",unlist(Tmax_postN_pct[10,]), pct="99")
  dev.off()
  # 
  # #Tmean
  pdf(file=paste0("Tmean_",file))
  # 
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmean_preN,"pred.Tmean_preN",unlist(Tmean_preN_pct[10,]), pct="99")
  # 
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[2,]), pct="1")
  lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[3,]), pct="5")
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[4,]), pct="10")
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[5,]), pct="20")
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[7,]), pct="80")
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[8,]), pct="90")
  lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[9,]), pct="95")
  # lapply(1:4,plot_pred_time,pred.Tmean_postN,"pred.Tmean_postN",unlist(Tmean_postN_pct[10,]), pct="99")
  dev.off()
  
}


# Effects at specific values ----
#matfit, matlow mathigh
# Tmax 4,7,8,9
# Tmin 5,10,11,12
# Tmean 6,13,14,15

#Low temperature: 2 3 4 5; 3 5th pct
#High temperature: 7 8 9 10; 9 95th pct


Tmean_mod_nb=c(6,13,14,15)
## Tmean_preN ----
{
  cumul_effects_Tmean_preN=function(mod_nb, 
                                    cen, at) {
    
    cp=crosspred(basis.Tmean_preN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmean_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmean_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  Tmean_preN_cumul_val=NULL
  ### Diff. temp ----
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmean_mod_nb[1],cumul_effects_Tmean_preN,cen=Tmean_preN_pct[6,],
                 at=unlist(Tmean_preN_pct[i,])))
    tmp=lapply(Tmean_mod_nb[1],cumul_effects_Tmean_preN,cen=Tmean_preN_pct[6,],
               at=unlist(Tmean_preN_pct[i,]))
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmean_preN_cumul_val=rbind(Tmean_preN_cumul_val,tmp)
    }
  }
  
  
}

## Tmean_postN ----
{
  cumul_effects_Tmean_postN=function(mod_nb, 
                                     cen, at) {
    
    cp=crosspred(basis.Tmean_postN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmean_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmean_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  Tmean_postN_cumul_val=NULL
  ### Diff. temp ----
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmean_mod_nb[1],cumul_effects_Tmean_postN,cen=Tmean_postN_pct[6,],
                 at=unlist(Tmean_postN_pct[i,])))
    tmp=lapply(Tmean_mod_nb[1],cumul_effects_Tmean_postN,cen=Tmean_postN_pct[6,],
               at=unlist(Tmean_postN_pct[i,]))
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmean_postN_cumul_val=rbind(Tmean_postN_cumul_val,tmp)
    }
  }
  
}


Tmax_mod_nb=c(4,7,8,9)
## Tmax_preN ----
{
  cumul_effects_Tmax_preN=function(mod_nb, 
                                   cen, at) {
    
    cp=crosspred(basis.Tmax_preN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmax_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmax_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  Tmax_preN_cumul_val=NULL
  ### Diff. temp ----
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmax_mod_nb[1],cumul_effects_Tmax_preN,cen=Tmax_preN_pct[6,],
                 at=unlist(Tmax_preN_pct[i,])))
    tmp=lapply(Tmax_mod_nb[1],cumul_effects_Tmax_preN,cen=Tmax_preN_pct[6,],
               at=unlist(Tmax_preN_pct[i,]))
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmax_preN_cumul_val=rbind(Tmax_preN_cumul_val,tmp)
    }
  }
  
  
}

## Tmax_postN ----
{
  cumul_effects_Tmax_postN=function(mod_nb, 
                                    cen, at) {
    
    cp=crosspred(basis.Tmax_postN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmax_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmax_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  Tmax_postN_cumul_val=NULL
  ### Diff. temp ----
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmax_mod_nb[1],cumul_effects_Tmax_postN,cen=Tmax_postN_pct[6,],
                 at=unlist(Tmax_postN_pct[i,])))
    tmp=lapply(Tmax_mod_nb[1],cumul_effects_Tmax_postN,cen=Tmax_postN_pct[6,],
               at=unlist(Tmax_postN_pct[i,]))
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmax_postN_cumul_val=rbind(Tmax_postN_cumul_val,tmp)
    }
  }
  
}


Tmin_mod_nb=c(5,10,11,12)
## Tmin_preN ----
{
  cumul_effects_Tmin_preN=function(mod_nb, 
                                   cen, at) {
    
    cp=crosspred(basis.Tmin_preN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmin_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmin_preN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=30:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  ### Diff. temp ----
  Tmin_preN_cumul_val=NULL
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmin_mod_nb[1],cumul_effects_Tmin_preN,cen=Tmin_preN_pct[6,],
                 at=unlist(Tmin_preN_pct[i,])))
    tmp=lapply(Tmin_mod_nb[1],cumul_effects_Tmin_preN,cen=Tmin_preN_pct[6,],
               at=unlist(Tmin_preN_pct[i,])) 
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmin_preN_cumul_val=rbind(Tmin_preN_cumul_val,tmp)
    }
  }
  
  
}

## Tmin_postN ----
{
  cumul_effects_Tmin_postN=function(mod_nb, 
                                    cen, at) {
    
    cp=crosspred(basis.Tmin_postN,mod[[mod_nb]], cen = cen, bylag=1, at = at)
    
    #where matlow and high have the same signs, find the min and max lags (x1, x2)
    where_min_max=data.frame(low=cp$matRRlow %>% as.vector, 
                             high=cp$matRRhigh %>% as.vector) %>%
      mutate(neg=case_when(low<1 & high < 1 ~"low")) %>% 
      mutate(pos=case_when(low>1 & high > 1 ~"high")) 
    min(which(where_min_max$neg=="low")) -> NEG_A
    max(which(where_min_max$neg=="low")) -> NEG_B
    
    min(which(where_min_max$pos=="high")) -> POS_A
    max(which(where_min_max$pos=="high")) -> POS_B
    
    #Nothing significant => exit
    NEG_CUM=NULL;POS_CUM=NULL
    which_lags_neg=NULL; which_lags_pos=NULL
    if(is.infinite(NEG_A) & is.infinite(POS_A)) return(NA)
    if(!is.infinite(NEG_A)) {
      cp=crosspred(basis.Tmin_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(NEG_A-1,NEG_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      NEG_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_neg=c(all_lags[NEG_B],all_lags[NEG_A])
      
    }
    
    if(!is.infinite(POS_A)) {
      cp=crosspred(basis.Tmin_postN, mod[[mod_nb]], cen = cen, at = at, lag=c(POS_A-1,POS_B-1))
      est=cp$allRRfit
      lower=cp$allRRlow
      upper=cp$allRRhigh
      POS_CUM=round(data.frame(est, lower, upper),3)
      
      #correspondance between lag and week (reverse)
      all_lags=91:1
      which_lags_pos=c(all_lags[POS_B],all_lags[POS_A])
      
    }
    
    return(list(which_lags_neg,which_lags_pos,NEG_CUM,POS_CUM))
  }
  
  Tmin_postN_cumul_val=NULL
  ### Diff. temp ----
  for (i in c(2:4,8:10)) {
    print(paste0("at: ",i))
    print(lapply(Tmin_mod_nb[1],cumul_effects_Tmin_postN,cen=Tmin_postN_pct[6,],
                 at=unlist(Tmin_postN_pct[i,])))
    tmp=lapply(Tmin_mod_nb[1],cumul_effects_Tmin_postN,cen=Tmin_postN_pct[6,],
               at=unlist(Tmin_postN_pct[i,]))
    if(is.list(do.call(c,tmp))) {
      tmp=cbind(i,do.call(cbind,do.call(c,do.call(c,tmp))))
      Tmin_postN_cumul_val=rbind(Tmin_postN_cumul_val,tmp)
    }
  }
  
}


Tmean_preN_cumul_val
Tmean_postN_cumul_val

Tmax_preN_cumul_val
Tmax_postN_cumul_val

Tmin_preN_cumul_val
Tmin_postN_cumul_val
