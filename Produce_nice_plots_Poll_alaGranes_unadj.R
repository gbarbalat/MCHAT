# hdr ----
#This script needs to be run after the main analysis script anal_data_NoStratif...

library(ggplot2)
library(ggplotify)
library(patchwork)

#1 is for unadjusted models, 4 adjusted on Tmean, 2 Tmax, 3 Tmin
which_Poll="PM10"

study="Overall"#Sex or Overall
if (which_Poll=="PM2.5") {
  Poll_mod_nb <- c(1)
} else if (which_Poll =="PM10") { 
  Poll_mod_nb <- c(2)
}else if (which_Poll =="NO2") {
  Poll_mod_nb <- c(3)}

preN_y1=0.99;preN_y2=1.015;#preN_y1=0.935;preN_y2=1.040;
postN_y1=0.98;postN_y2=1.010;


preN_xlim=c(29,0)
postN_xlim=c(90,0)
at_preN=c(29,seq(25,0,-5))#tickmarks
at_postN=c(90,seq(80,0,-10))#tickmarks
# at_preN1=29;at_preN2=0;
# at_postN1=90;at_postN2=0;

preN_xbreaks=c(1,seq(5,30,5))
postN_xbreaks=c(1,seq(10,91,10))

#### End of hdr
basis.Poll_preN=get(paste0("basis.",which_Poll,"_preN"))
basis.Poll_postN=get(paste0("basis.",which_Poll,"_postN"))

pred.Poll_preN=get(paste0("pred.",which_Poll,"_preN"))
pred.Poll_postN=get(paste0("pred.",which_Poll,"_postN"))

preN_breaks=c(preN_y1,(preN_y1+preN_y2)/2,preN_y2)
preN_ylim=c(preN_y1,preN_y2)

postN_breaks=c(postN_y1,(postN_y1+postN_y2)/2,postN_y2)
postN_ylim=c(postN_y1,postN_y2)


p=list()

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
                     ci.level)
  } else if (which_Poll=="PM10" & pre_or_post=="preN") {
    basis.PM10_preN=basis.Poll
    new_cp=crosspred(basis=basis.PM10_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag,
                     ci.level)
    
  } else if (which_Poll=="NO2" & pre_or_post=="preN") {
    basis.NO2_preN=basis.Poll
    new_cp=crosspred(basis=basis.NO2_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag,
                     ci.level)
    
  } else if (which_Poll=="PM2.5" & pre_or_post=="postN") {
    basis.PM2.5_postN=basis.Poll
    new_cp=crosspred(basis=basis.PM2.5_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag,
                     ci.level)
  } else if (which_Poll=="PM10" & pre_or_post=="postN") {
    basis.PM10_postN=basis.Poll
    new_cp=crosspred(basis=basis.PM10_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag,
                     ci.level)
    
  } else if (which_Poll=="NO2" & pre_or_post=="postN") {
    basis.NO2_postN=basis.Poll
    new_cp=crosspred(basis=basis.NO2_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag,
                     ci.level)
    
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

#### end of hdr

which_cp=1 #unadj


## preN ----
# Identify lags where CI doesn't include 1
significant_lags <- which(pred.Poll_preN[[which_cp]]$matRRlow["10",] > 1 | 
                            pred.Poll_preN[[which_cp]]$matRRhigh["10",] < 1)-1
#sort out their different ranges
ranges <- find_ranges(significant_lags)

#make it ggplot object for future patchwork
p[[1]] <- as.ggplot(function() {
 
   #main plot
   plot(pred.Poll_preN[[1]], 
       "slices", 
       var = 10,
       ylab = '', 
       xlab = '', 
       xlim=preN_xlim, xaxt='n',
       ylim=preN_ylim,yaxt='n',
       col = "black", 
       lwd = 2, 
       ci="line",
       ci.arg = list(col = "black")
  )
  axis(side = 1, at=at_preN,labels = preN_xbreaks)       
  axis(side = 2, preN_breaks) 

# compute for each set of range
new_cp=list()
if (!is.null(ranges)) {
  
for (i_range in 1:nrow(ranges)) {
  
  new_cp[[i_range]]=new_crosspred_Poll(
    basis.Poll=basis.Poll_preN,
    which_Poll = which_Poll, pre_or_post = "preN",
    model=mod[[Poll_mod_nb[which_cp]]],
    cen=0,
    at=10,
    lag=ranges[i_range,]
  )
  
  lines(new_cp[[i_range]], 
   "slices", 
   var = 10,
   ylab = '', 
   xlab = '', 
   xlim=preN_xlim, xaxt='n',
   ylim=preN_ylim,yaxt='n',
   col = "black", 
   lwd = 2,
   ci=ifelse(ranges[,1] == ranges[,2], "bar", "area"),
   ci.arg = list(col = "black")
  ) 
}
}
}
) %>% make_it_nice_gg_preN()+
  theme(plot.title = element_text(colour = "black"))


## `postN ----

# Identify lags where CI doesn't include 1
significant_lags <- which(pred.Poll_postN[[which_cp]]$matRRlow["10",] > 1 | 
                            pred.Poll_postN[[which_cp]]$matRRhigh["10",] < 1)-1
#sort out their different ranges
ranges <- find_ranges(significant_lags)

#make it ggplot object for future patchwork
p[[2]] <- as.ggplot(function() {
  
  #main plot
  plot(pred.Poll_postN[[which_cp]], 
       "slices", 
       var = 10,
       ylab = '', 
       xlab = '', 
       xlim=postN_xlim, xaxt='n',
       ylim=postN_ylim,yaxt='n',
       col = "black", 
       lwd = 2, 
       ci="line",
       ci.arg = list(col = "black")
  )
  axis(side = 1, at=at_postN,labels = postN_xbreaks)       
  axis(side = 2, postN_breaks) 
  
  # compute for each set of range
  new_cp=list()
  if (!is.null(ranges)) {
    
    for (i_range in 1:nrow(ranges)) {
      
      new_cp[[i_range]]=new_crosspred_Poll(
        basis.Poll=basis.Poll_postN,
        which_Poll = which_Poll, pre_or_post = "postN",
        model=mod[[Poll_mod_nb[which_cp]]],
        
        cen=0,
        at=10,
        lag=ranges[i_range,]
      )
      
      lines(new_cp[[i_range]], 
            "slices", 
            var = 10,
            ylab = '', 
            xlab = '', 
            xlim=postN_xlim, xaxt='n',
            ylim=postN_ylim,yaxt='n',
            col = "black", 
            lwd = 2,
            ci=ifelse(ranges[,1] == ranges[,2], "bar", "area"),
            ci.arg = list(col = "black")
      ) 
    }
  }
}
) %>% make_it_nice_gg_postN()+
  theme(plot.title = element_text(colour = "black"))

#combine ----
final_p <- p[[1]]+p[[2]]+patchwork::plot_layout(
  axis_titles = "collect", guides = "collect",nrow=1,ncol=2) 

print(final_p)
# 
# p=patchwork::wrap_plots(p,ncol=2 );
# #p=p & ylim(0.93,1.07)
# p


