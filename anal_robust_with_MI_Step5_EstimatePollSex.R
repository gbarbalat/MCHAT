rm(list=ls())
library(dplyr)
library(dlnm)
library(mice)
library(tidyr)
library(ggplot2)
library(patchwork)

#load imputed dataset
load("Imputed_MCHAT.RData")
ci.level <- 0.975
magic_nb <- 0.0125; #0.025 or 0.0125
alpha <- 0.025

# define cx basis
#prenatal and postnatal lags
nlags_preN<-29; i_lag_max_preN <- 30; step_preN<-7
nlags_postN<-23; i_lag_max_postN <- 24; step_postN<-30
nlags_postN<-19; i_lag_max_postN <- 20; step_postN<-30 #Discard last 4 months of exposure
nlags_postN<-90; i_lag_max_postN <- 91; step_postN<-7 #weekly 1 year and 3/4

#df for cx basis
grid1_Temp<-2#dose-resp relationship for temperature (ns function with INTERCEPT = FALSE)
grid2_Temp<-3 #lag-response relationship (ns function with INTERCEPT = TRUE)
grid1_Poll=NULL#dose-resp relationship for pollution; df=NULL (linear function)
grid2_Poll_preN=3
grid2_Poll_postN=4

#mice parameter
m <- 10
maxit <- 10

#exclude SEXE_ENF from all_covar
all_covar <- c(
  
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

# Ref temp values from original script
load("ref_temperature.RData")


#beginning of big function
#produces one plot per temperature indicator
#for heat or cold
plot_imputed_dlnm <- function(focus_on, adj) {
  
  
  focus_on_preN <-paste0(focus_on, "_preN_")#"PM2.5_preN_" 
  focus_on_postN <- paste0(focus_on, "_postN_")#"PM2.5_postN_" 
  
  adj_preN <-paste0(adj, "_preN_")#"Tmin_preN_" 
  adj_postN <- paste0(adj, "_postN_")#"Tmin_postN_"
  
  #ggplot parameters
  if (adj == "Tmean") {
    colour_title <- "darkgrey"
  } else if (adj == "Tmax") {
    colour_title <- "red"
  } else if (adj == "Tmin") {
    colour_title <- "darkblue"
  }
  
  
  # function define cx basis
  define_cx_b_AIC <-function(df_for_cb,SEX,var,varfun,nlags,i_lag_max,step, 
                             grid1, grid2) {
    i_lag<-1
    i_d<-1
    xx<-list()
    while (i_lag<=i_lag_max) {
      x<-list()
      for (i_step in 1:step) {
        var_i_d<-ifelse(grepl("preN", var, fixed = TRUE),i_d-1+365,i_d-1)
        x[[i_step]]<-df_for_cb  %>%
          select(paste0(var,var_i_d))
        i_d<-i_d+1
      }
      #xx=do.call(cbind,x)
      xx[[i_lag]]<-do.call(cbind,x) %>% rowMeans()
      
      i_lag<-i_lag+1
    }
    xxx<-do.call(cbind,xx)
    head(xxx)
    colSums(is.na(xxx))
    
    if(varfun=="lin") {argvar<-list(fun = varfun)
    } else {argvar<-list(fun = varfun, df = grid1)}
    
    
    ###
    sub<-xxx[1,];exphist(sub, times=length(sub), lag=length(sub)-1)
    xxx<-t(apply(xxx,1, function(sub) exphist(sub, 
                                              times=length(sub),
                                              lag=length(sub)-1)
    )
    )
    xxx[1,]
    ###
    cx_b<-crossbasis(x = xxx , 
                     lag = nlags, 
                     argvar = argvar, #ER rel
                     arglag = list(fun = "ns", df = grid2) #LR rel
                     
    )
    #colSums(is.na(basis.Tmax))
  }
  
  
  # fn to run dlnm analysis for each imputed dataset
  fn_dlnm_imp <- function(i,SEX) {
    
    # define formula for modelling
    formula <- as.formula(paste0("as.integer(ALL) ~ 
                           basis.Temp_preN + basis.Temp_postN + basis.Poll_preN + basis.Poll_postN +",
                                 paste0(all_covar, collapse = "+")))
    
    
    print(i)
    
    merged_imputed_i <- merged_imputed %>% complete(i) %>% select(- contains(focus_on))
    merged_ignore_withExp$SEXE_ENF <- merged_imputed_i$SEXE_ENF
    merged_ignore_withExp_Sex <- merged_ignore_withExp %>% 
      filter(SEXE_ENF==SEX)
    data_imp <- merged_imputed %>% complete(i) %>% select( - contains(focus_on)) %>%
      filter(SEXE_ENF==SEX)
    
    print(dim(merged_ignore_withExp_Sex))
    print(dim(data_imp))
    print(colSums(is.na(data_imp)))
    
    basis.Temp_preN<-define_cx_b_AIC(merged_ignore_withExp_Sex,SEX,var=adj_preN, #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                     varfun="ns",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                     grid1=grid1_Temp,grid2=grid2_Temp
                                     
    )
    
    basis.Temp_postN<-define_cx_b_AIC(merged_ignore_withExp_Sex,SEX,var=adj_postN,
                                      varfun="ns",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                      grid1=grid1_Temp,grid2=grid2_Temp
                                      
    )
    
    basis.Poll_preN=define_cx_b_AIC(merged_ignore_withExp_Sex,SEX,var=focus_on_preN, # NO2_preN_ PM2.5_preN_ PM10_preN_
                                    varfun="lin",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                    grid1=grid1_Poll,grid2=grid2_Poll_preN
                                    
    )
    
    basis.Poll_postN=define_cx_b_AIC(merged_ignore_withExp_Sex,SEX,var=focus_on_postN, # NO2_postN PM2.5_postN PM10_postN
                                     varfun="lin",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                     grid1=grid1_Poll,grid2=grid2_Poll_postN
                                     
    )
    
    ## run model on imputed dataset
    model <- glm(family=quasipoisson(),#poisson gaussian quasipoisson
                 #formula[[i]],
                 formula,
                 na.action = na.omit,
                 data_imp
    )
    print(model)
    
    ## cross pred 
    pred.Poll_preN <-crosspred(basis.Poll_preN,model,cen=0,bylag = 1,cumul = TRUE, at=10,ci.level = ci.level)
    pred.Poll_postN <-crosspred(basis.Poll_postN,model,cen=0,bylag = 1,cumul = TRUE, at=10,ci.level = ci.level)
    
    return(list(pred.Poll_preN=pred.Poll_preN, pred.Poll_postN=pred.Poll_postN ))
    
  }
  
  moderation_stat_SE <- function(Q1, SE1, Q2, SE2) {
    
    est=Q1-Q2
    se=sqrt(SE1^2+SE2^2)
    
    return(list(est=est, se=se))
    
  }
  
  #Now perform interaction tests
  fn_pseudo_IA <- function(i, results_miint_Male, results_miint_Female){
    
    #males.pred.Tmean_preN[[which_cp]]$matRRfit
    males.pred.Poll_preN <- results_miint_Male[[i]]$pred.Poll_preN
    females.pred.Poll_preN <- results_miint_Female[[i]]$pred.Poll_preN
    
    #calculate between group diff at each lag
    diff.cp_preN <- moderation_stat_SE(Q1=males.pred.Poll_preN$matfit,
                                       SE1=males.pred.Poll_preN$matse,
                                       Q2=females.pred.Poll_preN$matfit,
                                       SE2 = females.pred.Poll_preN$matse
    )
    
    males.pred.Poll_postN <- results_miint_Male[[i]]$pred.Poll_postN
    females.pred.Poll_postN <- results_miint_Female[[i]]$pred.Poll_postN
    
    #calculate between group diff at each lag
    diff.cp_postN <- moderation_stat_SE(Q1=males.pred.Poll_postN$matfit,
                                        SE1=males.pred.Poll_postN$matse,
                                        Q2=females.pred.Poll_postN$matfit,
                                        SE2 = females.pred.Poll_postN$matse
    )
    
    # #where are the lags that do not include 1 in their CI
    # significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
    #                             diff.cp$matRRhigh["10",] < 1)-1
    # #sort out their different ranges
    # ranges <- find_ranges(significant_lags)
    
    return(list(diff.cp_preN=diff.cp_preN,diff.cp_postN=diff.cp_postN))
  }
  
  # pool the results using Rubin's rule for each gp
  fn_Rubin_by_group <- function(results_miint,pre_post,fit, se) {
    
    results_sorted <- list()
    results_sorted$fit <- do.call(rbind, lapply(results_miint, function(x) x[[pre_post]][[fit]]))
    # Bind all $se vectors by rows
    results_sorted$se <- do.call(rbind, lapply(results_miint, function(x) x[[pre_post]][[se]]))
    # results_sorted$fit/se is a matrix (rows = imputations, columns = estimates/se per lag)
    
    ave_est <- colMeans(results_sorted$fit)
    var_inter <- sapply(seq_len(ncol(results_sorted$fit)), function(col) {
      (1/(m-1)) * sum((results_sorted$fit[, col] - ave_est[col])^2)
    })
    var_intra <- colMeans(results_sorted$se^2)
    total_var <- var_intra + (1 + 1/m)*var_inter
    r<-(1+1/m)*var_inter/var_intra
    vvv<-(m-1)*(1+1/r)^2 
    tal<-qt(magic_nb,df=vvv,lower.tail=F) 
    lCI <- ave_est - tal*sqrt(total_var)
    uCI <- ave_est + tal*sqrt(total_var) 
    pval<-2*(1-pt(q=abs(ave_est/sqrt(total_var)),df=vvv))  
    res_MI <-data.frame(cbind(ave_est,total_var,lCI,uCI,pval)) 
  }
  
  # pull results using Rubin rules between gp
  fn_Rubin <- function(results_miint,pre_post, fit, se) {
    results_sorted <- list()
    results_sorted$fit <- do.call(rbind, lapply(results_miint, function(x) x[[pre_post]][[fit]]))
    # Bind all $se vectors by rows
    results_sorted$se <- do.call(rbind, lapply(results_miint, function(x) x[[pre_post]][[se]]))
    # results_sorted$fit/se is a matrix (rows = imputations, columns = estimates/se per lag)
    
    ave_est <- colMeans(results_sorted$fit)
    var_inter <- sapply(seq_len(ncol(results_sorted$fit)), function(col) {
      (1/(m-1)) * sum((results_sorted$fit[, col] - ave_est[col])^2)
    })
    var_intra <- colMeans(results_sorted$se^2)
    total_var <- var_intra + (1 + 1/m)*var_inter
    r<-(1+1/m)*var_inter/var_intra
    vvv<-(m-1)*(1+1/r)^2 
    tal<-qt(magic_nb,df=vvv,lower.tail=F) 
    lCI <- ave_est - tal*sqrt(total_var)
    uCI <- ave_est + tal*sqrt(total_var) 
    pval<-2*(1-pt(q=abs(ave_est/sqrt(total_var)),df=vvv))  
    res_MI <-data.frame(cbind(ave_est,total_var,lCI,uCI,pval)) 
  }
  
  results_miint_Male <- lapply(1:m,fn_dlnm_imp, SEX="Male")
  results_miint_Female <- lapply(1:m,fn_dlnm_imp, SEX="Female")
  
  
  #Rubin per gp----
  within_Male_preN <- fn_Rubin_by_group(results_miint=results_miint_Male,pre_post="pred.Poll_preN",fit="matfit",se="matse")
  within_Male_postN <- fn_Rubin_by_group(results_miint=results_miint_Male,pre_post="pred.Poll_postN",fit="matfit",se="matse")
  
  within_Female_preN <- fn_Rubin_by_group(results_miint=results_miint_Female,pre_post="pred.Poll_preN",fit="matfit",se="matse")
  within_Female_postN <- fn_Rubin_by_group(results_miint=results_miint_Female,pre_post="pred.Poll_postN",fit="matfit",se="matse")
  
  #interaction test ----
  results_miint <- lapply(1:m, fn_pseudo_IA, results_miint_Male, results_miint_Female)
  
  res_MI_preN <- fn_Rubin(results_miint, pre_post="diff.cp_preN",fit="est", se="se")
  res_MI_postN <- fn_Rubin(results_miint, pre_post="diff.cp_postN",fit="est", se="se")
  
  
  plot_anal_robust <- function(res_MI_Male, res_MI_Female, pre_post) {
    
    plot_title <- ""#paste0("Adj. for ", adj)
    
    
    if (pre_post == "preN") {
      x_breaks <- c(1, 5, 10, 15, 20, 25, 30)
    } else if (pre_post == "postN") {
      x_breaks <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
    }
    
    x_lab <- ifelse(pre_post=="preN","Weeks since conception (Prenatal period)",
                    "Weeks since delivery (Postnatal period)" )
    
    #combine the data
    res_MI_Male$Sex <- "Male"
    res_MI_Male$Lag <-rev(seq_len(nrow(res_MI_Male)))
    
    
    res_MI_Female$Sex <- "Female"
    res_MI_Female$Lag <- rev(seq_len(nrow(res_MI_Male)))
    
    
    res_MI <- rbind(res_MI_Male, res_MI_Female)
    
    #get results of interaction test lag by lag and make it the pval column
    res_MI_IA <- get(paste0("res_MI_",pre_post))
    res_MI$pval <- res_MI_IA$pval
    
    res_MI <- res_MI[rev(seq_len(nrow(res_MI))), ]
    res_MI <- res_MI |>
      dplyr::mutate(
        RR = exp(ave_est),
        RR_lCI = exp(lCI),
        RR_uCI = exp(uCI)
      ) %>%
      mutate(significant = pval <= alpha)
    
    res_MI <- res_MI %>%
      arrange(Sex, Lag) %>%
      group_by(Sex) %>%
      mutate(
        sig_group = {
          is_sig <- significant
          rle_sig <- rle(is_sig)
          group_ids <- rep(NA_integer_, length(is_sig))
          group_counter <- 1
          idx <- 1
          for (i in seq_along(rle_sig$lengths)) {
            len <- rle_sig$lengths[i]
            val <- rle_sig$values[i]
            if (val) {
              group_ids[idx:(idx + len - 1)] <- group_counter
              group_counter <- group_counter + 1
            }
            idx <- idx + len
          }
          group_ids
        }
      ) %>%
      ungroup()
    
    ggplot(res_MI, aes(x = Lag, y = RR, color = Sex, fill = Sex)) +
      
      scale_x_continuous(breaks = x_breaks) +
      
      # Shaded ribbon only where pval <= 0.05
      # geom_ribbon(
      #   data = res_MI %>% filter(significant),
      #   aes(ymin = RR_lCI, ymax = RR_uCI,fill = Sex),
      #    alpha = 0.3,  color = NA
      #   ) +
      
      geom_ribbon(
        data = res_MI %>% filter(significant),
        aes(ymin = RR_lCI, ymax = RR_uCI, fill = Sex, group = interaction(Sex, sig_group)),
        alpha = 0.3,
        color = NA
      ) +
      
      
      # Main line
      geom_line( ) +
      #geom_point() +
      
      # Dashed lines for lCI and uCI
      geom_line(aes(y = RR_lCI), linetype = "dashed") +#, color = colour_signif
      geom_line(aes(y = RR_uCI), linetype = "dashed") +#, color = colour_signif
      
      # Reference line
      geom_hline(yintercept = 1, color = "black") +
      
      # Labels and styling
      labs(x = x_lab, y = "Relative Risk", title = plot_title) +
      
      # Custom colors for male/female
      scale_color_manual(values = c("Male" = "darkgreen", "Female" = "purple")) +
      scale_fill_manual(values = c("Male" = "darkgreen", "Female" = "purple")) +
      
      theme(
        legend.position = "none",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 13, angle = 90, face="bold"),
        axis.text.y = element_text(size = 11, angle = 90),
        axis.title.x = element_text(size = 13, face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(face = "bold", size = 15, 
                                  hjust = 0.5, vjust = -5),
        axis.line = element_line(color = "black")
      )
    
  }
  
  preN_Poll <- plot_anal_robust(res_MI_Male=within_Male_preN, 
                                res_MI_Female=within_Female_preN, 
                                pre_post="preN")
  
  postN_Poll <- plot_anal_robust(res_MI_Male=within_Male_postN, 
                                 res_MI_Female=within_Female_postN, 
                                 pre_post="postN")
  
  return(list(preN_Poll=preN_Poll, postN_Poll=postN_Poll))
  
  #end of big function
}

#plots for PM2.5
PM2.5_Tmean_plots <- plot_imputed_dlnm(focus_on="PM2.5", adj="Tmean")
PM2.5_Tmax_plots <- plot_imputed_dlnm(focus_on="PM2.5", adj="Tmax")
PM2.5_Tmin_plots <- plot_imputed_dlnm(focus_on="PM2.5", adj="Tmin")

Poll_Tmean_preN  <- PM2.5_Tmean_plots$preN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
PM2.5_Tmean_preN  <- PM2.5_Tmean_plots$preN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmean_postN <- PM2.5_Tmean_plots$postN_Poll + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
PM2.5_Tmean_postN <- PM2.5_Tmean_plots$postN_Poll + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmax_preN  <- PM2.5_Tmax_plots$preN_Poll  + theme(axis.title.x = element_blank())
Poll_Tmax_postN  <- PM2.5_Tmax_plots$postN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmin_preN   <- PM2.5_Tmin_plots$preN_Poll   + theme(axis.title.y = element_blank())
Poll_Tmin_postN  <- PM2.5_Tmin_plots$postN_Poll  + theme(axis.title.y = element_blank())

(Poll_Tmean_preN + Poll_Tmean_postN) /
  (Poll_Tmax_preN + Poll_Tmax_postN) /
  (Poll_Tmin_preN + Poll_Tmin_postN)


#plots for PM10
PM10_Tmean_plots <- plot_imputed_dlnm(focus_on="PM10", adj="Tmean")
PM10_Tmax_plots <- plot_imputed_dlnm(focus_on="PM10", adj="Tmax")
PM10_Tmin_plots <- plot_imputed_dlnm(focus_on="PM10", adj="Tmin")

Poll_Tmean_preN  <- PM10_Tmean_plots$preN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
PM10_Tmean_preN  <- PM10_Tmean_plots$preN_Poll  + theme(axis.title.x = element_blank())
Poll_Tmean_postN <- PM10_Tmean_plots$postN_Poll + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
PM10_Tmean_postN <- PM10_Tmean_plots$postN_Poll + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmax_preN  <- PM10_Tmax_plots$preN_Poll  + theme(axis.title.x = element_blank())
Poll_Tmax_postN  <- PM10_Tmax_plots$postN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmin_preN   <- PM10_Tmin_plots$preN_Poll   + theme(axis.title.y = element_blank())
Poll_Tmin_postN  <- PM10_Tmin_plots$postN_Poll  + theme(axis.title.y = element_blank())

(Poll_Tmean_preN + Poll_Tmean_postN) /
  (Poll_Tmax_preN + Poll_Tmax_postN) /
  (Poll_Tmin_preN + Poll_Tmin_postN)


#Plots for NO2
NO2_Tmean_plots <- plot_imputed_dlnm(focus_on="NO2", adj="Tmean")
NO2_Tmax_plots <- plot_imputed_dlnm(focus_on="NO2", adj="Tmax")
NO2_Tmin_plots <- plot_imputed_dlnm(focus_on="NO2", adj="Tmin")

Poll_Tmean_preN  <- NO2_Tmean_plots$preN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
NO2_Tmean_preN  <- NO2_Tmean_plots$preN_Poll  + theme(axis.title.y = element_blank())
Poll_Tmean_postN <- NO2_Tmean_plots$postN_Poll + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
NO2_Tmean_postN <- NO2_Tmean_plots$postN_Poll + theme(axis.title.y = element_blank())
Poll_Tmax_preN  <- NO2_Tmax_plots$preN_Poll  + theme(axis.title.x = element_blank())
Poll_Tmax_postN  <- NO2_Tmax_plots$postN_Poll  + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
Poll_Tmin_preN   <- NO2_Tmin_plots$preN_Poll   + theme(axis.title.y = element_blank())
Poll_Tmin_postN  <- NO2_Tmin_plots$postN_Poll  + theme(axis.title.y = element_blank())

(Poll_Tmean_preN + Poll_Tmean_postN) /
  (Poll_Tmax_preN + Poll_Tmax_postN) /
  (Poll_Tmin_preN + Poll_Tmin_postN)

#for another article format ...
(PM2.5_Tmean_preN + PM2.5_Tmean_postN) /
  (PM10_Tmean_preN + PM10_Tmean_postN) /
  (NO2_Tmean_preN + NO2_Tmean_postN)
