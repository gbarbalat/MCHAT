
# function to calculate btw group differences
moderation_stat_ratio=function(Q1, lower1, upper1, Q2, lower2, upper2) {
  # SE1=(upper1-lower1)/3.92
  # SE2=(upper2-lower2)/3.92
  SElog1=(log(upper1)-log(lower1))/3.92
  SElog2=(log(upper2)-log(lower2))/3.92
  est=Q1/Q2
  lower=exp(log(Q1/Q2)-1.96*sqrt(SElog1^2+SElog2^2))
  upper=exp(log(Q1/Q2)+1.96*sqrt(SElog1^2+SElog2^2))
  
  out=data.frame(est=est, lower=lower, upper=upper)
  
  return(out)
}

# The following sees if there are any signif differences over the whole study period (prenat OR postnat)
moderation_stat=function(Q1, lower1, upper1, Q2, lower2, upper2) {
  SElog1=(log(upper1)-log(lower1))/3.92
  SElog2=(log(upper2)-log(lower2))/3.92
  est=data.frame(round(Q1/Q2,3))
  lower=exp(log(Q1/Q2)-1.96*sqrt(SElog1^2+SElog2^2))
  upper=exp(log(Q1/Q2)+1.96*sqrt(SElog1^2+SElog2^2))
  
  return(list(matRRfit=est, matRRlow=lower,matRRhigh= upper))
}

calculate_btwgpdiff <- function(Exp_idx, which_temp, pre_post) {
  
  lag_to_return=NULL
  
  # Models concerned by Tmean
  Tmean_mod_nb=c(6,13,14,15)
  # Models concerned by Tmin
  Tmin_mod_nb=c(5,10,11,12)
  # Models concerned by Tmax
  Tmax_mod_nb=c(4,7,8,9)
  
  #cross-pred at stake (here always the first one from the model unadj for Poll)
  which_cp=1
  
  # from lags in the dlnm sense to weeks 
  all_lags_preN=30:1
  all_lags_postN=91:1
  
  ## Tmean ----

  if (which_temp=="Tmean" & pre_post=="_preN") {
    
  print("Tmean_preN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmean_preN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmean_preN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmean_preN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmean_preN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmean_preN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmean_preN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmean_preN=function(basis.Tmean_preN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmean_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmean_preN(
        basis.Tmean_preN=males.basis.Tmean_preN,
        model=males.mod[[Tmean_mod_nb[which_cp]]],
        cen=Tmean_preN_pct[6,],
        at=Tmean_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmean_preN(
        basis.Tmean_preN=females.basis.Tmean_preN,
        model=females.mod[[Tmean_mod_nb[which_cp]]],
        cen=Tmean_preN_pct[6,],
        at=Tmean_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmean_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmean_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmean_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmean_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmean_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmean_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
        lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
      
    }
  }
  }
  
  
  
  ### Tmean_postN
  else if (which_temp=="Tmean" & pre_post=="_postN") {

  print("Tmean_postN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmean_postN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmean_postN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmean_postN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmean_postN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmean_postN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmean_postN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmean_postN=function(basis.Tmean_postN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmean_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmean_postN(
        basis.Tmean_postN=males.basis.Tmean_postN,
        model=males.mod[[Tmean_mod_nb[which_cp]]],
        cen=Tmean_postN_pct[6,],
        at=Tmean_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmean_postN(
        basis.Tmean_postN=females.basis.Tmean_postN,
        model=females.mod[[Tmean_mod_nb[which_cp]]],
        cen=Tmean_postN_pct[6,],
        at=Tmean_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmean_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmean_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmean_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmean_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmean_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmean_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
        lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
    }
  }
  
  }

  
  
  
  
  
  
  
  ## Tmax ----
  else if (which_temp=="Tmax" & pre_post=="_preN") {
  
  ### Tmax_preN
  print("Tmax_preN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmax_preN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmax_preN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmax_preN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmax_preN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmax_preN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmax_preN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmax_preN=function(basis.Tmax_preN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmax_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmax_preN(
        basis.Tmax_preN=males.basis.Tmax_preN,
        model=males.mod[[Tmax_mod_nb[which_cp]]],
        cen=Tmax_preN_pct[6,],
        at=Tmax_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmax_preN(
        basis.Tmax_preN=females.basis.Tmax_preN,
        model=females.mod[[Tmax_mod_nb[which_cp]]],
        cen=Tmax_preN_pct[6,],
        at=Tmax_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmax_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmax_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmax_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmax_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmax_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmax_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
        lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
    }
  }
  }
  
  
  ### Tmax_postN
  else if (which_temp=="Tmax" & pre_post=="_postN") {
  
  print("Tmax_postN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmax_postN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmax_postN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmax_postN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmax_postN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmax_postN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmax_postN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmax_postN=function(basis.Tmax_postN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmax_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmax_postN(
        basis.Tmax_postN=males.basis.Tmax_postN,
        model=males.mod[[Tmax_mod_nb[which_cp]]],
        cen=Tmax_postN_pct[6,],
        at=Tmax_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmax_postN(
        basis.Tmax_postN=females.basis.Tmax_postN,
        model=females.mod[[Tmax_mod_nb[which_cp]]],
        cen=Tmax_postN_pct[6,],
        at=Tmax_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmax_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmax_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmax_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmax_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmax_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmax_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
        lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
    }
  }
  
  }
  
  ## Tmin ----

  
  ### Tmin_preN
  else if (which_temp=="Tmin" & pre_post=="_preN") {
    
  print("Tmin_preN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmin_preN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmin_preN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmin_preN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmin_preN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmin_preN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmin_preN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmin_preN=function(basis.Tmin_preN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmin_preN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmin_preN(
        basis.Tmin_preN=males.basis.Tmin_preN,
        model=males.mod[[Tmin_mod_nb[which_cp]]],
        cen=Tmin_preN_pct[6,],
        at=Tmin_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmin_preN(
        basis.Tmin_preN=females.basis.Tmin_preN,
        model=females.mod[[Tmin_mod_nb[which_cp]]],
        cen=Tmin_preN_pct[6,],
        at=Tmin_preN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmin_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmin_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmin_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmin_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmin_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmin_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
        lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
    }
  }
  }
  
  ### Tmin_postN
  else if (which_temp=="Tmin" & pre_post=="_postN") {

  print("Tmin_postN")
  #calculate between group diff at each lag
  diff.cp=moderation_stat(Q1=males.pred.Tmin_postN[[which_cp]]$matRRfit,
                          lower1=males.pred.Tmin_postN[[which_cp]]$matRRlow,
                          upper1 = males.pred.Tmin_postN[[which_cp]]$matRRhigh,
                          Q2=females.pred.Tmin_postN[[which_cp]]$matRRfit,
                          lower2=females.pred.Tmin_postN[[which_cp]]$matRRlow,
                          upper2 = females.pred.Tmin_postN[[which_cp]]$matRRhigh
  )
  
  #where are the lags that do not include 1 in their CI
  significant_lags <- which(diff.cp$matRRlow[Exp_idx,] > 1 | 
                              diff.cp$matRRhigh[Exp_idx,] < 1)-1
  #sort out their different ranges
  ranges <- find_ranges(significant_lags)
  
  # compute new crosspred where significant differences for each lag
  new_crosspred_Tmin_postN=function(basis.Tmin_postN,model,cen,at,lag) {
    
    new_cp=crosspred(basis=basis.Tmin_postN,
                     model=model,
                     cen=cen,
                     at=at,
                     lag = lag)
    
  }
  
  # compute within- and between-group diff for each set of lag range
  new_cp_males=list()
  new_cp_females=list()
  moder=list()
  if (!is.null(ranges)) {
    for (i_range in 1:nrow(ranges)) {
      
      #Where is the lag?
      print("lag")
      print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
      
      
      #within males
      new_cp_males[[i_range]]=new_crosspred_Tmin_postN(
        basis.Tmin_postN=males.basis.Tmin_postN,
        model=males.mod[[Tmin_mod_nb[which_cp]]],
        cen=Tmin_postN_pct[6,],
        at=Tmin_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("males")
      print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                       lower=round(new_cp_males[[i_range]]$allRRlow,3),
                       upper=round(new_cp_males[[i_range]]$allRRhigh,3))
      )
      
      #within females
      new_cp_females[[i_range]]=new_crosspred_Tmin_postN(
        basis.Tmin_postN=females.basis.Tmin_postN,
        model=females.mod[[Tmin_mod_nb[which_cp]]],
        cen=Tmin_postN_pct[6,],
        at=Tmin_postN_pct[Exp_idx,] %>% unlist,
        lag=ranges[i_range,]
      )
      print("females")
      print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                       lower=round(new_cp_females[[i_range]]$allRRlow,3),
                       upper=round(new_cp_females[[i_range]]$allRRhigh,3))
      )
      
      
      #between groups
      moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.Tmin_preN.matRRfit,
                                             lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.Tmin_preN.matRRlow,
                                             upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.Tmin_preN.matRRhigh,
                                             Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.Tmin_preN.matRRfit,
                                             lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.Tmin_preN.matRRlow,
                                             upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.Tmin_preN.matRRhigh
      )
      print("Btw")
      print(moder[[i_range]])
      
      if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
         lag_to_return=rbind(lag_to_return,ranges[i_range,])
        
      }
      
    }
  }
  }
  
  return(lag_to_return)
  
}




calculate_btwgpdiff_Poll <- function(which_Poll, pre_post,which_cp) {
  
  lag_to_return=NULL
  
  # Models concerned by PM2.5
  PM2.5_mod_nb=c(1,8,11,14)
  # Models concerned by PM10
  PM10_mod_nb=c(2,9,12,15)
  # Models concerned by NO2
  NO2_mod_nb=c(3,7,10,13)
  
  # from lags in the dlnm sense to weeks 
  all_lags_preN=30:1
  all_lags_postN=91:1
    
    ## PM2.5 ----
    
    if (which_Poll=="PM2.5" & pre_post=="_preN") {
      
    print("PM2.5_preN")
    
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,]
    )
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_PM2.5_preN=function(basis.PM2.5_preN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.PM2.5_preN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_PM2.5_preN(
          basis.PM2.5_preN=males.basis.Poll_preN,
          model=males.mod[[PM2.5_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_PM2.5_preN(
          basis.PM2.5_preN=females.basis.Poll_preN,
          model=females.mod[[PM2.5_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.PM2.5_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.PM2.5_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.PM2.5_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.PM2.5_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.PM2.5_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.PM2.5_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    
    else if (which_Poll=="PM2.5" & pre_post=="_postN") {
      
    ### PM2.5_postN
    print("PM2.5_postN")
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,]
    )
    
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_PM2.5_postN=function(basis.PM2.5_postN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.PM2.5_postN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_PM2.5_postN(
          basis.PM2.5_postN=males.basis.Poll_postN,
          model=males.mod[[PM2.5_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_PM2.5_postN(
          basis.PM2.5_postN=females.basis.Poll_postN,
          model=females.mod[[PM2.5_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.PM2.5_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.PM2.5_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.PM2.5_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.PM2.5_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.PM2.5_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.PM2.5_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    
    ## PM10 ----
    
    else if (which_Poll=="PM10" & pre_post=="_preN") {
      
    print("PM10_preN")
    
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,]
    )
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_PM10_preN=function(basis.PM10_preN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.PM10_preN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_PM10_preN(
          basis.PM10_preN=males.basis.Poll_preN,
          model=males.mod[[PM10_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_PM10_preN(
          basis.PM10_preN=females.basis.Poll_preN,
          model=females.mod[[PM10_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.PM10_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.PM10_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.PM10_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.PM10_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.PM10_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.PM10_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    
    else if (which_Poll=="PM10" & pre_post=="_postN") {
      
    ### PM10_postN
    print("PM10_postN")
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,]
    )
    
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_PM10_postN=function(basis.PM10_postN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.PM10_postN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_PM10_postN(
          basis.PM10_postN=males.basis.Poll_postN,
          model=males.mod[[PM10_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_PM10_postN(
          basis.PM10_postN=females.basis.Poll_postN,
          model=females.mod[[PM10_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.PM10_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.PM10_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.PM10_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.PM10_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.PM10_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.PM10_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    
    ## NO2 ----
  
    else if (which_Poll=="NO2" & pre_post=="_preN") {
      
    print("NO2_preN")
    
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_preN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_preN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_preN[[which_cp]]$matRRhigh[1:10,]
    )
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_NO2_preN=function(basis.NO2_preN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.NO2_preN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_preN[ranges[i_range,1]+1],all_lags_preN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_NO2_preN(
          basis.NO2_preN=males.basis.Poll_preN,
          model=males.mod[[NO2_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_NO2_preN(
          basis.NO2_preN=females.basis.Poll_preN,
          model=females.mod[[NO2_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.NO2_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.NO2_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.NO2_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.NO2_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.NO2_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.NO2_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    else if (which_Poll=="NO2" & pre_post=="_postN") {
      
    ### NO2_postN
    print("NO2_postN")
    #calculate between group diff at each lag
    diff.cp=moderation_stat(Q1=males.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower1=males.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper1 =males.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,],
                            Q2=females.pred.Poll_postN[[which_cp]]$matRRfit[1:10,],
                            lower2=females.pred.Poll_postN[[which_cp]]$matRRlow[1:10,],
                            upper2 = females.pred.Poll_postN[[which_cp]]$matRRhigh[1:10,]
    )
    
    
    #where are the lags that do not include 1 in their CI
    significant_lags <- which(diff.cp$matRRlow["10",] > 1 | 
                                diff.cp$matRRhigh["10",] < 1)-1
    #sort out their different ranges
    ranges <- find_ranges(significant_lags)
    
    # compute new crosspred where significant differences for each lag
    new_crosspred_NO2_postN=function(basis.NO2_postN,model,cen,at,lag) {
      
      new_cp=crosspred(basis=basis.NO2_postN,
                       model=model,
                       cen=cen,
                       at=at,
                       lag = lag)
      
    }
    
    # compute within- and between-group diff for each set of lag range
    new_cp_males=list()
    new_cp_females=list()
    moder=list()
    if (!is.null(ranges)) {
      for (i_range in 1:nrow(ranges)) {
        
        #Where is the lag?
        print("lag")
        print(c(all_lags_postN[ranges[i_range,1]+1],all_lags_postN[ranges[i_range,2]+1]) %>% rev)
        
        
        #within males
        new_cp_males[[i_range]]=new_crosspred_NO2_postN(
          basis.NO2_postN=males.basis.Poll_postN,
          model=males.mod[[NO2_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("males")
        print(data.frame(est=round(new_cp_males[[i_range]]$allRRfit,3),
                         lower=round(new_cp_males[[i_range]]$allRRlow,3),
                         upper=round(new_cp_males[[i_range]]$allRRhigh,3))
        )
        
        #within females
        new_cp_females[[i_range]]=new_crosspred_NO2_postN(
          basis.NO2_postN=females.basis.Poll_postN,
          model=females.mod[[NO2_mod_nb[which_cp]]],
          cen=0,
          at=10,
          lag=ranges[i_range,]
        )
        print("females")
        print(data.frame(est=round(new_cp_females[[i_range]]$allRRfit,3),
                         lower=round(new_cp_females[[i_range]]$allRRlow,3),
                         upper=round(new_cp_females[[i_range]]$allRRhigh,3))
        )
        
        
        #between groups
        moder[[i_range]]=moderation_stat_ratio(Q1=new_cp_males[[i_range]]$allRRfit,#males.pred.NO2_preN[[which_cp]]$matRRfit,
                                               lower1=new_cp_males[[i_range]]$allRRlow,#males.pred.NO2_preN[[which_cp]]$matRRlow,
                                               upper1 =new_cp_males[[i_range]]$allRRhigh, #males.pred.NO2_preN[[which_cp]]$matRRhigh,
                                               Q2=new_cp_females[[i_range]]$allRRfit,#females.pred.NO2_preN[[which_cp]]$matRRfit,
                                               lower2=new_cp_females[[i_range]]$allRRlow,#females.pred.NO2_preN[[which_cp]]$matRRlow,
                                               upper2 = new_cp_females[[i_range]]$allRRhigh#females.pred.NO2_preN[[which_cp]]$matRRhigh
        )
        print("Btw")
        print(moder[[i_range]])
        if (moder[[i_range]]$lower > 1 | moder[[i_range]]$upper < 1) {
          lag_to_return=rbind(lag_to_return,ranges[i_range,])
          
        }
        
      }
    }
    }
    
    return(lag_to_return)
  }
  