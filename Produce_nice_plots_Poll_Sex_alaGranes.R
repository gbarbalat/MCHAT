# hdr ----
rm(list=ls())

library(dlnm)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggplotify)
library(patchwork)

path <- "./37WA/" #""  "./37WA/"

source("calculate_btwgpdiff.R")
add_on <- "" #"_season" ""
which_Poll="PM10"#PM2.5 PM10 NO2
ci.level <- 0.975; magic_nb <- 2.24; #1.96 or 2.2414

#which_cp=; varies 2 3 4
#1 is for unadjusted models, 3 adjusted on PM2.5, 4 PM10, 2 NO2
# 9 is for p95 and 3 is for p5
study="Sex"#Sex or Overall


preN_y1=0.975;preN_y2=1.025;
#no2
if (which_Poll=="NO2") {
postN_y1=0.985;postN_y2=1.010;
}
#pm
if (which_Poll=="PM2.5" | which_Poll=="PM10") {
  postN_y1=0.965;postN_y2=1.030;
}

preN_breaks=c(preN_y1,(preN_y1+preN_y2)/2,preN_y2)
preN_ylim=c(preN_y1,preN_y2)
postN_breaks=c(postN_y1,(postN_y1+postN_y2)/2,postN_y2)
postN_ylim=c(postN_y1,postN_y2)

preN_breaks=c(preN_y1,(preN_y1+preN_y2)/2,preN_y2)
preN_ylim=c(preN_y1,preN_y2)

postN_breaks=c(postN_y1,(postN_y1+postN_y2)/2,postN_y2)
postN_ylim=c(postN_y1,postN_y2)


#code title 
title_tmin <- "" # "Adj. for Tmin"
title_tmax <- "" # "Adj. for Tmax"
title_tmean <- "" # "Adj. for Tmean"


# # Models concerned by Tmean
# Tmean_mod_nb=c(6,13,14,15)
# # Models concerned by Tmin
# Tmin_mod_nb=c(5,10,11,12)
# # Models concerned by Tmax
# Tmax_mod_nb=c(4,7,8,9)

# formula=list(formula_PM2.5,formula_PM10,formula_NO2,
#              formula_Tmax,formula_Tmin,formula_Tmean,
#              formula_Tmax_NO2,formula_Tmax_PM2.5,formula_Tmax_PM10,
#              formula_Tmin_NO2,formula_Tmin_PM2.5,formula_Tmin_PM10,
#              formula_Tmean_NO2,formula_Tmean_PM2.5,formula_Tmean_PM10)

if (which_Poll=="PM2.5") {
  Poll_mod_nb <- c(1,8,11,14)
  } else if (which_Poll =="PM10") { 
  Poll_mod_nb <- c(2,9,12,15)
  }else if (which_Poll =="NO2") {
  Poll_mod_nb <- c(3,7,10,13)}

#### End of hdr


#load and transform ----
load(paste0(path,"Male",add_on,".RData"))

males.pred.Poll_preN=get(paste0("pred.",which_Poll,"_preN"))
males.pred.Poll_postN=get(paste0("pred.",which_Poll,"_postN"))
males.basis.Poll_preN=get(paste0("basis.",which_Poll,"_preN"))
males.basis.Poll_postN=get(paste0("basis.",which_Poll,"_postN"))
males.mod=mod

load(paste0(path,"Female",add_on,".RData"))

females.pred.Poll_preN=get(paste0("pred.",which_Poll,"_preN"))
females.pred.Poll_postN=get(paste0("pred.",which_Poll,"_postN"))
females.basis.Poll_preN=get(paste0("basis.",which_Poll,"_preN"))
females.basis.Poll_postN=get(paste0("basis.",which_Poll,"_postN"))
females.mod=mod

preN_xlim=c(nlags_preN,0)
postN_xlim=c(nlags_postN,0)
at_preN=c(nlags_preN,seq(i_lag_max_preN-5,0,-5))#tickmarks
at_postN=c(nlags_postN,seq(nlags_postN-10,0,-10))#tickmarks

preN_xbreaks=c(1,seq(5,i_lag_max_preN,5))
postN_xbreaks=c(1,seq(10,i_lag_max_postN,10))

# function to compute new crosspred
new_crosspred_Poll=function(basis.Poll, which_Poll,pre_or_post,
                            model,cen,at,lag, ci.level) {
  

  if (which_Poll=="PM2.5" & pre_or_post=="preN") {
    basis.PM2.5_preN=basis.Poll
    new_cp=crosspred(basis=basis.PM2.5_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
  } else if (which_Poll=="PM10" & pre_or_post=="preN") {
    basis.PM10_preN=basis.Poll
    new_cp=crosspred(basis=basis.PM10_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
    
  } else if (which_Poll=="NO2" & pre_or_post=="preN") {
    basis.NO2_preN=basis.Poll
    new_cp=crosspred(basis=basis.NO2_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
    
  } else if (which_Poll=="PM2.5" & pre_or_post=="postN") {
    basis.PM2.5_postN=basis.Poll
    new_cp=crosspred(basis=basis.PM2.5_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
  } else if (which_Poll=="PM10" & pre_or_post=="postN") {
    basis.PM10_postN=basis.Poll
    new_cp=crosspred(basis=basis.PM10_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
    
  } else if (which_Poll=="NO2" & pre_or_post=="postN") {
    basis.NO2_postN=basis.Poll
    new_cp=crosspred(basis=basis.NO2_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag, 
                     ci.level=ci.level)
    
  }
}



#function to find range of significant lags
find_ranges <- function(x) {
  if (length(x) == 0) return(NULL)
  
  gaps <- diff(x) > 1
  start <- c(x[1], x[which(gaps) + 1])
  end <- c(x[which(gaps)], x[length(x)])
  
  return(cbind(start, end))
}

make_it_nice_gg_preN = function(initial_plot) {
  initial_plot + labs(x="Weeks since conception (Prenatal period)", y="Relative Risk") +
    theme(legend.position = "none",
          legend.text=element_text(size=11),
          axis.text = element_blank(),
          axis.title.y = element_text(size = 15, angle = 90, vjust = NULL),#element_blank(),#element_text(size = 15),
          axis.title.x = element_text(size = 15,vjust=NULL),#element_blank(),#element_text(size = 15),
          strip.text = element_text(size=13),
          panel.background = element_rect(fill = "white",colour = "white"),
          strip.background = element_rect(fill = "white",colour = "black"),
          plot.title = element_text(colour = "black",face="bold", size=15, 
                                    hjust=0.5,vjust=-5) ,
          #aspect.ratio = 0.5
    )
}

make_it_nice_gg_postN = function(initial_plot) {
  initial_plot + labs(x="Weeks since delivery (Postnatal period)", y="Relative Risk") +
    theme(legend.position = "none",
          legend.text=element_text(size=11),
          axis.text = element_blank(),
          axis.title.y = element_text(size = 15, angle = 90, vjust = NULL),#element_blank(),#element_text(size = 15),
          axis.title.x = element_text(size = 15,vjust=NULL),#element_blank(),#element_text(size = 15),
          strip.text = element_text(size=13),
          panel.background = element_rect(fill = "white",colour = "white"),
          strip.background = element_rect(fill = "white",colour = "black"),
          plot.title = element_text(colour = "black",face="bold", size=15, 
                                    hjust=0.5,vjust=-5) ,
          #aspect.ratio = 0.5
    )
}

# # function to find significant lags
# signif_lags_preN <- function(which_cp) {
#   #calculate between group diff at each lag
#   diff.cp=moderation_stat(Q1=males.pred.Poll_preN[[which_cp]]$matRRfit["10",],
#                           lower1=males.pred.Poll_preN[[which_cp]]$matRRlow["10",],
#                           upper1 = males.pred.Poll_preN[[which_cp]]$matRRhigh["10",],
#                           Q2=females.pred.Poll_preN[[which_cp]]$matRRfit["10",],
#                           lower2=females.pred.Poll_preN[[which_cp]]$matRRlow["10",],
#                           upper2 = females.pred.Poll_preN[[which_cp]]$matRRhigh["10",]
#   )
#   
#   #where are the lags that do not include 1 in their CI
#   significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
#                               diff.cp$matRRhigh["10",] < 1)-1
#   
#   ranges <- find_ranges(significant_lags)
# }

#### end of hdr


#prepare list of graphs
p=list()


# Tmin ----
which_cp=3

## preN ----

#main plot
p[[5]] <- as.ggplot(function() {
    plot.crosspred(males.pred.Poll_preN[[which_cp]], #pred.Poll_preN[[3]], 
                   "slices", 
                   var = 10,
                   ylab = '', 
                   xlab = '', 
                   xlim=preN_xlim, xaxt='n',
                   ylim=preN_ylim,yaxt='n',
                   col = "darkgreen", 
                   lwd = 2, lty=1,
                   ci="line",
                   ci.level = ci.level,
                   
                   ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
    )
    lines(females.pred.Poll_preN[[which_cp]],  
          "slices", 
          var = 10,
          # ylab = '', 
          # xlab = '', 
          # xlim=preN_xlim, xaxt='n',
          # ylim=preN_ylim,yaxt='n',
          col = "purple", 
          lwd = 2, lty=1,
          ci="line",
          ci.level = ci.level,
          
          ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
    )
    axis(side = 1, at=at_preN,labels = preN_xbreaks)       
    axis(side = 2, preN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_preN",which_cp=which_cp, ci.level = ci.level) 
  

    # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_preN()+
    ggtitle(paste0(title_tmin)
            
  ) +
  theme(plot.title = element_text(colour = "blue"))


## postN ----
p[[6]] <- as.ggplot(function() {
  #main
  plot.crosspred(males.pred.Poll_postN[[which_cp]], #pred.Poll_preN[[3]], 
                 "slices", 
                 var = 10,
                 ylab = '', 
                 xlab = '', 
                 xlim=postN_xlim, xaxt='n',
                 ylim=postN_ylim,yaxt='n',
                 col = "darkgreen", 
                 lwd = 2, lty=1,
                 ci="line",
                 ci.level = ci.level,
                 
                 ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
  )
  lines(females.pred.Poll_postN[[which_cp]],  
        "slices", 
        var = 10,
        col = "purple", 
        lwd = 2, lty=1,
        ci="line",
        ci.level = ci.level,
        
        ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
  )
  axis(side = 1, at=at_postN,labels = postN_xbreaks)       
  axis(side = 2, postN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_postN",which_cp=which_cp, ci.level = ci.level) 
  
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_postN()+
  ggtitle(paste0(title_tmin)
          
  ) +
  theme(plot.title = element_text(colour = "blue"))



# Tmax ----
which_cp=2

## preN ----

#main plot
p[[3]] <- as.ggplot(function() {
  plot.crosspred(males.pred.Poll_preN[[which_cp]], #pred.Poll_preN[[3]], 
                 "slices", 
                 var = 10,
                 ylab = '', 
                 xlab = '', 
                 xlim=preN_xlim, xaxt='n',
                 ylim=preN_ylim,yaxt='n',
                 col = "darkgreen", 
                 lwd = 2, lty=1,
                 ci="line",
                 ci.level = ci.level,
                 
                 ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
  )
  lines(females.pred.Poll_preN[[which_cp]],  
        "slices", 
        var = 10,
        # ylab = '', 
        # xlab = '', 
        # xlim=preN_xlim, xaxt='n',
        # ylim=preN_ylim,yaxt='n',
        col = "purple", 
        lwd = 2, lty=1,
        ci="line",
        ci.level = ci.level,
        
        ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
  )
  axis(side = 1, at=at_preN,labels = preN_xbreaks)       
  axis(side = 2, preN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_preN",which_cp=which_cp, ci.level = ci.level) 
  
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_preN()+
  ggtitle(paste0(title_tmax)
          
  ) +
  theme(plot.title = element_text(colour = "red"))


## postN ----
p[[4]] <- as.ggplot(function() {
  #main
  plot.crosspred(males.pred.Poll_postN[[which_cp]], #pred.Poll_preN[[3]], 
                 "slices", 
                 var = 10,
                 ylab = '', 
                 xlab = '', 
                 xlim=postN_xlim, xaxt='n',
                 ylim=postN_ylim,yaxt='n',
                 col = "darkgreen", 
                 lwd = 2, lty=1,
                 ci="line",
                 ci.level = ci.level,
                 
                 ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
  )
  lines(females.pred.Poll_postN[[which_cp]],  
        "slices", 
        var = 10,
        col = "purple", 
        lwd = 2, lty=1,
        ci="line",
        ci.level = ci.level,
        
        ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
  )
  axis(side = 1, at=at_postN,labels = postN_xbreaks)       
  axis(side = 2, postN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_postN",which_cp=which_cp, ci.level = ci.level) 
  
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_postN()+
  ggtitle(paste0(title_tmax)
          
  ) +
  theme(plot.title = element_text(colour = "red"))





# Tmean ----
which_cp=4

## preN ----

#main plot
p[[1]] <- as.ggplot(function() {
  plot.crosspred(males.pred.Poll_preN[[which_cp]], #pred.Poll_preN[[3]], 
                 "slices", 
                 var = 10,
                 ylab = '', 
                 xlab = '', 
                 xlim=preN_xlim, xaxt='n',
                 ylim=preN_ylim,yaxt='n',
                 col = "darkgreen", 
                 lwd = 2, lty=1,
                 ci="line",
                 ci.level = ci.level,
                 
                 ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
  )
  lines(females.pred.Poll_preN[[which_cp]],  
        "slices", 
        var = 10,
        # ylab = '', 
        # xlab = '', 
        # xlim=preN_xlim, xaxt='n',
        # ylim=preN_ylim,yaxt='n',
        col = "purple", 
        lwd = 2, lty=1,
        ci="line",
        ci.level = ci.level,
        
        ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
  )
  axis(side = 1, at=at_preN,labels = preN_xbreaks)       
  axis(side = 2, preN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_preN",which_cp=which_cp, ci.level = ci.level) 
  
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_preN,
        which_Poll = which_Poll, pre_or_post = "preN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_preN()+
  ggtitle(paste0(title_tmean)
          
  ) +
  theme(plot.title = element_text(colour = "darkgrey"))


## postN ----
p[[2]] <- as.ggplot(function() {
  #main
  plot.crosspred(males.pred.Poll_postN[[which_cp]], #pred.Poll_preN[[3]], 
                 "slices", 
                 var = 10,
                 ylab = '', 
                 xlab = '', 
                 xlim=postN_xlim, xaxt='n',
                 ylim=postN_ylim,yaxt='n',
                 col = "darkgreen", 
                 lwd = 2, lty=1,
                 ci="line",
                 ci.level = ci.level,
                 
                 ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
  )
  lines(females.pred.Poll_postN[[which_cp]],  
        "slices", 
        var = 10,
        col = "purple", 
        lwd = 2, lty=1,
        ci="line",
        ci.level = ci.level,
        
        ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
  )
  axis(side = 1, at=at_postN,labels = postN_xbreaks)       
  axis(side = 2, postN_breaks)       
  
  lag_returned <- calculate_btwgpdiff_Poll(which_Poll=which_Poll, pre_post="_postN",which_cp=which_cp, ci.level = ci.level) 
  
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(lag_returned)) {
    for (i_range in 1:nrow(lag_returned)) {
      # Loop code here
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=males.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=males.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "darkgreen", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("darkgreen", alpha.f = 0.3))
      )
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=females.basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=females.mod[[Poll_mod_nb[which_cp]]],
        cen=0,
        at=10,
        lag=lag_returned[i_range,], 
        ci.level=ci.level
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            # ylab = '', 
            # xlab = '', 
            # xlim=preN_xlim, xaxt='n',
            # ylim=preN_ylim,yaxt='n',
            col = "purple", 
            lwd = 2, lty=1,
            ci="area",
            ci.level = ci.level,
            
            ci.arg = list(col = adjustcolor("purple", alpha.f = 0.3))
      )
      
      
    }
  } 
}
) %>% make_it_nice_gg_postN()+
  ggtitle(paste0(title_tmean)
          
  ) +
  theme(plot.title = element_text(colour = "darkgrey"))




#combine ----
final_p <- p[[1]]+p[[2]]+
  #p[[3]]+p[[4]]+
  #p[[5]]+p[[6]]+
  patchwork::plot_layout(
  axis_titles = "collect", guides = "collect",
  #nrow=3,ncol=2
  nrow=1,ncol=2
  ) 

print(final_p)
# 
# p=patchwork::wrap_plots(p,ncol=2 );
# #p=p & ylim(0.93,1.07)
# p

# check values for adjustment of axes lims
summary(males.pred.Poll_preN[[4]]$matRRlow["10",])
summary(males.pred.Poll_preN[[4]]$matRRhigh["10",])
summary(females.pred.Poll_preN[[4]]$matRRlow["10",])
summary(females.pred.Poll_preN[[4]]$matRRhigh["10",])

summary(males.pred.Poll_postN[[4]]$matRRlow["10",])
summary(males.pred.Poll_postN[[4]]$matRRhigh["10",])
summary(females.pred.Poll_postN[[4]]$matRRlow["10",])
summary(females.pred.Poll_postN[[4]]$matRRhigh["10",])
