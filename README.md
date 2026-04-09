# MCHAT
This the code for the MCHAT study

# Paper ref

# Matching scripts and analyses  
Examples are given for the pollution analyses but same principles need to be applied for the temperature analyses.  
Beware: ci.level=0.975 instead of 0.95 (see methods from the paper)

## Pollution unadjusted
1- anal_data_NoStratif_RefTempAsIan_unadj.R  
2- keep obtained crosspred (e.g. pred.PM2.5_postN)  
3- Produce_nice_plots_Poll_alaGranes_unadj.R  
4- Copy paste R plots onto odg file, rearrange, save as and import into your paper

## Pollution adjusted   
1- anal_data_NoStratif_RefTempAsIan.R  
2- keep obtained crosspred (e.g. pred.PM2.5_postN)  
3- Produce_nice_plots_Poll_alaGranes.R  

## Pollution E/R  
1- anal_data_NoStratif_RefTempAsIan_ERcurve.R  

## Pollution Sex-stratified  
1- anal_data_StratifSex_RefTempAsIan.R, save Male.RData and Female.RData  
2- Produce_nice_plots_Poll_Sex_alagranes.R

## Pollution adjusted when increasing prenatal exposure window

## Pollution Sex-stratified when increasing prenatal exposure window

## Pollution adjusted MI

## Pollution Sex-stratified MI
