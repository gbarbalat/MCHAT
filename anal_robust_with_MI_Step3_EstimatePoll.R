rm(list=ls())
library(dplyr)
library(dlnm)
library(mice)
library(tidyr)
library(ggplot2)
library(patchwork)


#beginning of big function
#produces one plot per pollutant
#adj for each of the three temperature indicators
ci.level <- 0.975
magic_nb <- 0.0125; #0.025 or 0.0125
alpha <- 0.025

plot_imputed_dlnm <- function(focus_on, adj) {
  
  focus_on_preN <-paste0(focus_on, "_preN_")#"PM2.5_preN_" 
  focus_on_postN <- paste0(focus_on, "_postN_")#"PM2.5_postN_" 
  
  adj_preN <-paste0(adj, "_preN_")#"Tmin_preN_" 
  adj_postN <- paste0(adj, "_postN_")#"Tmin_postN_" 
  
  #ggplot parameters
  if (adj == "Tmean") {
    colour_signif <- "darkgrey"
    colour_NS <- "lightgrey"
  } else if (adj == "Tmax") {
    colour_signif <- "red"
    colour_NS <- "pink"
  } else if (adj == "Tmin") {
    colour_signif <- "darkblue"
    colour_NS <- "lightblue"
  }
  
  #mice parameter
  m <- 10
  maxit <- 10
  
  #load imputed dataset
  load("Imputed_MCHAT.RData")
  
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
  
  
  # Ref temp values
  #This is the ref values used in cross preds: median and various percentiles
  #one for each period (prenat vs. postnat) and each temp variable (Tmax, Tmin, Tmean)
  #from original script
  load("ref_temperature.RData")
  
  # function define cx basis
  define_cx_b_AIC <-function(var,varfun,nlags,i_lag_max,step, 
                             grid1, grid2) {
    i_lag<-1
    i_d<-1
    xx<-list()
    while (i_lag<=i_lag_max) {
      x<-list()
      for (i_step in 1:step) {
        var_i_d<-ifelse(grepl("preN", var, fixed = TRUE),i_d-1+365,i_d-1)
        x[[i_step]]<-merged_ignore_withExp %>% 
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
  
  
  Temp_preN_pct<-get(paste0(adj_preN, "pct"))
  Temp_postN_pct<-get(paste0(adj_postN, "pct"))
  
  basis.Temp_preN<-define_cx_b_AIC(var=adj_preN, #Tmean_preN_ Tmax_preN_ Tmin_preN_
                                   varfun="ns",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                   grid1=grid1_Temp,grid2=grid2_Temp
                                   
  )
  
  basis.Temp_postN<-define_cx_b_AIC(var=adj_postN,
                                    varfun="ns",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                    grid1=grid1_Temp,grid2=grid2_Temp
                                    
  )
  
  basis.Poll_preN=define_cx_b_AIC(var=focus_on_preN, # NO2_preN_ PM2.5_preN_ PM10_preN_
                                  varfun="lin",nlags=nlags_preN, i_lag_max = i_lag_max_preN, step=step_preN, 
                                  grid1=grid1_Poll,grid2=grid2_Poll_preN
                                  
  )
  
  basis.Poll_postN=define_cx_b_AIC(var=focus_on_postN, # NO2_postN PM2.5_postN PM10_postN
                                   varfun="lin",nlags= nlags_postN,i_lag_max = i_lag_max_postN, step=step_postN, 
                                   grid1=grid1_Poll,grid2=grid2_Poll_postN
                                   
  )
  
  # define formula for modelling
  formula <- as.formula(paste0("as.integer(ALL) ~ 
                           basis.Temp_preN + basis.Temp_postN + basis.Poll_preN + basis.Poll_postN +",
                               paste0(all_covar, collapse = "+")))
  
  # fn to run dlnm analysis for each imputed dataset
  fn_dlnm_imp <- function(i) {
    ## imputed dataset
    data_imp <- merged_imputed %>% complete(i) %>% select(-contains(adj), - contains(focus_on))# %>%  cbind(merged_[sort(keep_id$case),] %>% select(starts_with("Tmin")))
    
    ## run model on imputed dataset
    model <- glm(family=quasipoisson(),#poisson gaussian quasipoisson
                 #formula[[i]],
                 formula,
                 na.action = na.omit,
                 data_imp
    )
    
    ## cross pred 
    # pred.Temp_preN <-crosspred(basis.Temp_preN,model,cen=Temp_preN_pct[6,],bylag = 1,cumul = TRUE, at=Temp_preN_pct %>% unlist)
    # pred.Temp_postN <-crosspred(basis.Temp_postN,model,cen=Temp_postN_pct[6,],bylag = 1,cumul = TRUE, at=Temp_postN_pct %>% unlist)
    
    pred.Poll_preN <-crosspred(basis.Poll_preN,model,cen=0,bylag = 1,cumul = TRUE, at=10, ci.level=ci.level)
    pred.Poll_postN <-crosspred(basis.Poll_postN,model,cen=0,bylag = 1,cumul = TRUE, at=10, ci.level=ci.level)
    
    RRfit_preN <- pred.Poll_preN$matRRfit
    RRlow_preN <- pred.Poll_preN$matRRlow
    RRhigh_preN <- pred.Poll_preN$matRRhigh
    
    fit_preN <- pred.Poll_preN$matfit
    se_preN <- pred.Poll_preN$matse
    
    RRfit_postN <- pred.Poll_postN$matRRfit
    RRlow_postN <- pred.Poll_postN$matRRlow
    RRhigh_postN <- pred.Poll_postN$matRRhigh
    
    fit_postN <- pred.Poll_postN$matfit
    se_postN <- pred.Poll_postN$matse
    
    return(list(fit_preN=fit_preN, se_preN=se_preN, fit_postN=fit_postN, se_postN=se_postN ))
    
  }
  results_miint <- lapply(1:m,fn_dlnm_imp)
  
  # pool the results using Rubin's rule 
  
  # fn Bind all $fit vectors by rows
  fn_Rubin <- function(fit, se) {
    results_sorted <- list()
    results_sorted$fit <- do.call(rbind, lapply(results_miint, function(x) x[[fit]]))
    # Bind all $se vectors by rows
    results_sorted$se <- do.call(rbind, lapply(results_miint, function(x) x[[se]]))
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
  
  res_MI_preN <- fn_Rubin(fit="fit_preN", se="se_preN")
  res_MI_postN <- fn_Rubin(fit="fit_postN", se="se_postN")
  
  print(res_MI_preN)
  print(res_MI_postN)
  
  plot_anal_robust <- function(res_MI, pre_post) {
    
    
    plot_title <- ""#paste0("Adj. for ", adj)
    
    if (pre_post == "preN") {
      x_breaks <- c(1, 5, 10, 15, 20, 25, 30)
    } else if (pre_post == "postN") {
      x_breaks <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
    }
    
    x_lab <- ifelse(pre_post=="preN","Weeks since conception (Prenatal period)",
                    "Weeks since delivery (Postnatal period)" )
    
    res_MI <- res_MI[rev(seq_len(nrow(res_MI))), ]
    res_MI <- res_MI |>
      dplyr::mutate(
        RR = exp(ave_est),
        RR_lCI = exp(lCI),
        RR_uCI = exp(uCI),
        Lag = seq_len(nrow(res_MI))  # create lag variable if needed
      ) %>%
      mutate(significant = pval <= alpha)
    
    ggplot(res_MI, aes(x = Lag, y = RR)) +
      
      scale_x_continuous(breaks = x_breaks) +
      
      # Shaded ribbon only where pval <= 0.05
      geom_ribbon(
        data = res_MI %>% filter(significant),
        aes(ymin = RR_lCI, ymax = RR_uCI),
        fill = colour_NS
      ) +
      
      # Main line
      geom_line( colour = colour_signif) +
      #geom_point() +
      
      # Dashed lines for lCI and uCI
      # Dashed lines for lCI and uCI
      geom_line(aes(y = RR_lCI), linetype = "dashed", color = colour_signif) +
      geom_line(aes(y = RR_uCI), linetype = "dashed", color = colour_signif) +
      
      # Reference line
      geom_hline(yintercept = 1, color = "black") +
      
      # Labels and styling
      labs(x = x_lab, y = "Relative Risk", title = plot_title) +
      theme(
        legend.position = "none",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 13, angle = 90, face="bold"),
        axis.text.y = element_text(size = 11, angle = 90),
        axis.title.x = element_text(size = 13, face="bold"),
        panel.background = element_blank(),
        plot.title = element_text(colour = colour_signif, face = "bold", size = 15, 
                                  hjust = 0.5, vjust = -5),
        axis.line = element_line(color = "black")
      )
    
  }
  
  preN_Poll <- plot_anal_robust(res_MI_preN, "preN")
  postN_Poll <- plot_anal_robust(res_MI_postN, "postN")
  
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
