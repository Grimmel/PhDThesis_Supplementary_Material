library(terra)
library(kernlab)
library(dismo)
library(maxnet)
library(randomForest)
library(scales)


# Wilson, P. D. (2011). Distance-based methods for the analysis of maps produced by 
#     species distribution models. Methods in Ecology and Evolution, 2(6), 623–633. 
#     https://doi.org/10.1111/j.2041-210X.2011.00115.x
HellingerDist <- function (mat1,mat2)
{
  p1 <- sum(mat1,na.rm=T)
  p2 <- sum(mat2,na.rm=T)
  return(sqrt(0.5*sum((sqrt(mat1/p1) - sqrt(mat2/p2))^2,na.rm=T)))
}
getStack <- function(landscape){
  s <- rast(c(paste("F:\\PhD\\prediction_extent\\data\\ls\\species2\\LS",landscape,"_temp.asc",sep=""),
         paste("F:\\PhD\\prediction_extent\\data\\ls\\species2\\LS",landscape,"_precip.asc",sep=""),
         paste("F:\\PhD\\prediction_extent\\data\\ls\\species2\\LS",landscape,"_habitat.asc",sep="")))
  names(s) <-c('temp','precip','hab')
  s<- ifel(is.na(s), -99, s)
  return(s)
}
getSample <- function(n,stack,occ,t){
  if (t=='background'){
    backgroundpt <- spatSample(occ,10000,as.points=TRUE)
    backgroundpt$pa <- 0
    background <- terra::extract(stack,backgroundpt)
    return(cbind(backgroundpt,background))
  } else {
    samplept <- spatSample(occurrences,n,method='random',as.points=TRUE)
    names(samplept) <- c('pa')
    # Select occurrences only
    samplept <- samplept[samplept$pa==1,]
    sample <- terra::extract(stack,samplept)
    return(cbind(samplept,sample))
  }
}
genON <- function(presences,background,calStack,xferStack,calLS,xferLS,rep){
  ###### Create models for occupied niche and available environment #####
  #OCSVM
  ocsvm <- ksvm(pa~.,data=values(presences[,c('temp','precip','hab','pa')]),
                type="one-svc",kernel="rbfdot")
  ocsvm_e <- ksvm(pa~.,data=values(background[,c('temp','precip','hab','pa')]),
                  type="one-svc",kernel="rbfdot")
  oc_C <- predict(calStack,ocsvm,type="decision")
  oc_X <- predict(xferStack,ocsvm,type="decision")
  oc_BC <- predict(calStack,ocsvm_e,type="decision")
  oc_BX <- predict(xferStack,ocsvm_e,type="decision")
  # Isolation forest
  iforest <- isolation.forest(values(presences[,c('temp','precip','hab')]))
  iforest_e <- isolation.forest(values(background[,c('temp','precip','hab')]))
  if_C<-(predict(calStack,iforest,type="score"))
  if_X <-(predict(xferStack,iforest,type="score"))
  if_BX<-(predict(xferStack,iforest_e,type="score"))
  if_BC <-(predict(calStack,iforest_e,type="score"))
  # MESS (need to convert from terra rast to Raster stack)
  cStack <- stack(raster(calStack$temp),raster(calStack$precip),raster(calStack$hab))
  xStack <- stack(raster(xferStack$temp),raster(xferStack$precip),raster(xferStack$hab))
  ms_X <- mess(xStack,values(presences[,c('temp','precip','hab')]))
  ms_C <- mess(cStack,values(presences[,c('temp','precip','hab')]))
  ms_BX <- mess(xStack,values(background[,c('temp','precip','hab')]))
  ms_BC <- mess(cStack,values(background[,c('temp','precip','hab')]))
  
  ##### Convert continuous maps to binary outlier maps #####
  # Get map values within the calibration environment at sample locations
  ocSamp <- extract(oc_C,presences)
  ocESamp <- extract(oc_BC,background)
  ifSamp <- extract(if_C,presences)
  ifESamp <- extract(if_BC,background)
  msSamp <- extract(rast(ms_C),presences)
  msESamp <- extract(rast(ms_BC),background)
  # Obtain the threshold value for each map
  ocTH <-  min(ocSamp$lyr1)
  oceTH <- min(ocESamp$lyr1)
  ifTH <-  quantile(ifSamp$lyr1,0.95)
  ifeTH <- quantile(ifESamp$lyr1,0.95)
  # Convert to binary maps
  ocNi <- ifel(oc_X>ocTH,1,0)
  ocENi <- ifel(oc_BX>oceTH,1,0)
  ifNi <- ifel(if_X<ifTH,1,0)
  ifENi <- ifel(if_BX<ifeTH,1,0)
  msNi <- ifel(rast(ms_X)>0,1,0)
  crs(msNi) <- crs(ocNi)
  msENi <- ifel(rast(ms_BX)>0,1,0)
  crs(msENi) <- crs(ocENi)
  
  out <- rast(list(ocNi,ifNi,msNi,ocENi,ifENi,msENi))
  names(out) <- c('ocON','ifON','msON','ocEN','ifEN','msEN')
  writeRaster(ocNi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\oc_ON_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(ifNi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\if_ON_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(msNi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\ms_ON_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(ocENi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\oc_EN_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(ifENi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\if_EN_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(msENi,paste("F:\\PhD\\prediction_extent\\model_outputs\\profile\\ms_EN_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  return(out)
}
runMaxnetModel <- function(pa,calStack,xferStack,nicheStack,calLS,xferLS,rep){
  m1 <- maxnet(pa$pa,values(pa[,c('temp','precip','hab')]))

  e <- evaluate(values(pa[pa$pa==1,c('temp','precip','hab')]),values(pa[pa$pa==0,c('temp','precip','hab')]),m1,type="logistic",clamp=F)
  th <- threshold(e)
  xpred <- predict(xferStack,m1,type="logistic",clamp=F)
  npred <- predict(nicheStack,m1,type="logistic",clamp=F)
  binxpred <- ifel(xpred>th$spec_sens,1,0)
  out <- rast(list(xpred,binxpred))
  names(out) <- c("pred_suit","pred_occ")
  writeRaster(npred,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\max_niche_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(xpred,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\max_suit_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(binxpred,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\max_bin_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  return(out)
}
runRFModel <- function(pa,calStack,xferStack,nicheStack,calLS,xferLS,rep){
  pa <- as.data.frame(pa)
  # Random forest performs better with less background. Reduced to 1000 to match occurrence sampling intensity.
  bg <- pa[(pa$pa==0),]
  bg <- bg[sample(nrow(bg),1000),]
  paRF <- rbind(pa[(pa$pa==1),],bg)
  m <- randomForest(as.factor(pa)~temp+precip+hab,data=paRF)
  e <- evaluate(paRF[paRF$pa==1,c('temp','precip','hab')],paRF[paRF$pa==0,c('temp','precip','hab')],m,type="prob")
  th <- threshold(e)
  xpred <- predict(xferStack,m,type="prob")
  npred <- predict(nicheStack,m,type="prob")
  binxpred <- ifel(xpred$X1>th$spec_sens,1,0)
  out <- rast(list(xpred$X1,binxpred))
  names(out) <- c("pred_suit","pred_occ")
  writeRaster(npred$X1,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\rf_niche_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(xpred$X1,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\rf_suit_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  writeRaster(binxpred,paste("F:\\PhD\\prediction_extent\\model_outputs\\predictions3\\rf_bin_cLS",calLS,"_xLS",xferLS,"_",rep,".tif",sep=""))
  return(out)
}
calculateMetrics <- function(predictions,ocniche,suit,alg,rep,calLS,xferLS){
  analyse <-as.data.frame(rast(list(predictions,ocniche,suit)))
  analyse$calocc <- analyse$ocON+analyse$ifON+analyse$msON
  analyse$calen <- analyse$ocEN+analyse$ifEN+analyse$msEN
  analyse[analyse$calocc>0,]$calocc <- 1
  analyse[analyse$calen>0,]$calen <- 1
  
  analyse$pred_suit_scale <- scales::rescale(analyse$pred_suit, to=c(min(analyse$suitability),max(analyse$suitability)))
  analyse$pred_occ_scale <- scales::rescale(analyse$pred_suit, to=c(min(analyse$prob_occ),max(analyse$prob_occ)))
  analyse$suit_diff <- analyse$suitability-analyse$pred_suit_scale
  analyse$pocc_diff <- analyse$prob_occ-analyse$pred_occ_scale
  
  # Occupied niche - in/out
  # Environment - in/out
  # out ON, in env
  h_suit <- HellingerDist(analyse$pred_suit,analyse$suitability)
  hON_suit <- HellingerDist(analyse[analyse$calocc==1,]$pred_suit,analyse[analyse$calocc==1,]$suitability)
  hNON_suit <- HellingerDist(analyse[analyse$calocc==0,]$pred_suit,analyse[analyse$calocc==0,]$suitability)
  hENV_suit <- HellingerDist(analyse[analyse$calen==1,]$pred_suit,analyse[analyse$calen==1,]$suitability)
  hNENV_suit <- HellingerDist(analyse[analyse$calen==0,]$pred_suit,analyse[analyse$calen==0,]$suitability)
  hENV_NON_suit <- HellingerDist(analyse[(analyse$calen==1) & (analyse$calocc==0),]$pred_suit,analyse[(analyse$calen==1) & (analyse$calocc==0),]$suitability)
  
  h_poc <- HellingerDist(analyse$pred_suit,analyse$prob_occ)
  hON_poc <- HellingerDist(analyse[analyse$calocc==1,]$pred_suit,analyse[analyse$calocc==1,]$prob_occ)
  hNON_poc <- HellingerDist(analyse[analyse$calocc==0,]$pred_suit,analyse[analyse$calocc==0,]$prob_occ)
  hENV_poc <- HellingerDist(analyse[analyse$calen==1,]$pred_suit,analyse[analyse$calen==1,]$prob_occ)
  hNENV_poc <- HellingerDist(analyse[analyse$calen==0,]$pred_suit,analyse[analyse$calen==0,]$prob_occ)
  hENV_NON_poc <- HellingerDist(analyse[(analyse$calen==1) & (analyse$calocc==0),]$pred_suit,analyse[(analyse$calen==1) & (analyse$calocc==0),]$prob_occ)
  
  # RMSE
  rmseAll_suit <- sqrt(mean((analyse$suit_diff)^2))
  rmseON_suit <- sqrt(mean((analyse[analyse$calocc==1,]$suit_diff)^2))
  rmseNON_suit <- sqrt(mean((analyse[analyse$calocc==0,]$suit_diff)^2))
  rmseENV_suit <- sqrt(mean((analyse[analyse$calen==1,]$suit_diff)^2))
  rmseNENV_suit <- sqrt(mean((analyse[analyse$calen==0,]$suit_diff)^2))
  rmseENV_NON_suit <- sqrt(mean((analyse[(analyse$calen==1) & (analyse$calocc==0),]$suit_diff)^2))
  
  rmseAll_poc <- sqrt(mean((analyse$pocc_diff)^2))
  rmseON_poc<- sqrt(mean((analyse[analyse$calocc==1,]$pocc_diff)^2))
  rmseNON_poc <- sqrt(mean((analyse[analyse$calocc==0,]$pocc_diff)^2))
  rmseENV_poc <- sqrt(mean((analyse[analyse$calen==1,]$pocc_diff)^2))
  rmseNENV_poc <- sqrt(mean((analyse[analyse$calen==0,]$pocc_diff)^2))
  rmseENV_NON_poc <- sqrt(mean((analyse[(analyse$calen==1) & (analyse$calocc==0),]$pocc_diff)^2))
  
  # Occurrences
  predArea <- sum(analyse$pred_occ)
  predONArea <- sum(analyse[analyse$calocc==1,]$pred_occ)
  predNONArea <- sum(analyse[analyse$calocc==0,]$pred_occ)
  predENVArea <- sum(analyse[analyse$calen==1,]$pred_occ)
  predNENVArea <- sum(analyse[analyse$calen==0,]$pred_occ)
  predENV_NONArea <- sum(analyse[(analyse$calen==1) & (analyse$calocc==0),]$pred_occ)
  
  out <- data.frame(alg=alg,rep=rep,calLS=calLS,xferLS=xferLS,h_suit=h_suit,hON_suit=hON_suit,hNON_suit = hNON_suit,hENV_suit = hENV_suit,hNENV_suit = hNENV_suit,hENV_NON_suit = hENV_NON_suit,h_poc = h_poc,hON_poc = hON_poc,hNON_poc = hNON_poc,hENV_poc=hENV_poc,hNENV_poc = hNENV_poc,hENV_NON_poc = hENV_NON_poc,
           rmseAll_suit=rmseAll_suit,rmseON_suit=rmseON_suit,rmseNON_suit=rmseNON_suit,rmseENV_suit=rmseENV_suit,rmseNENV_suit=rmseNENV_suit,rmseENV_NON_suit=rmseENV_NON_suit,rmseAll_poc=rmseAll_poc,rmseON_poc=rmseON_poc,rmseNON_poc=rmseNON_poc,rmseENV_poc=rmseENV_poc,rmseNENV_poc=rmseNENV_poc,rmseENV_NON_poc=rmseENV_NON_poc,
           predArea=predArea,predONArea=predONArea,predNONArea=predNONArea,predENVArea=predENVArea,predNENVArea=predNENVArea,predENV_NONArea=predENV_NONArea)
  return(out)
}

lsPairs <- list(c(6,5),c(6,2),c(1,2),c(11,12),c(10,3),c(3,1),c(13,3),c(5,6),c(2,6),c(2,1),c(12,11),c(3,10),c(1,3),c(3,13))
nicheStack <- rast(c("F:\\PhD\\prediction_extent\\data\\test_temp.tif",
                     "F:\\PhD\\prediction_extent\\data\\test_precip.tif",
                     "F:\\PhD\\prediction_extent\\data\\test_hab.tif"))
names(nicheStack) <- c('temp','precip','hab')
results <- 0
for (pair in lsPairs){
  print(pair)
  calStack <- getStack(pair[1])
  xferStack <- getStack(pair[2])
  test <- rast(c(paste("F:\\PhD\\prediction_extent\\data\\ls\\species2\\LS",pair[2],"_suit.asc",sep= ""),
                 paste("F:\\PhD\\prediction_extent\\sims\\species2\\LS",pair[2],"_average.asc",sep="")))
  names(test) <- c('suitability','prob_occ')
  repCount <- 1
  for(sim in c(10,110,210,310,410,50,150,250,350,450)){
    occurrences <- rast(paste("F:\\PhD\\prediction_extent\\sims\\species2\\LS",pair[1],"_rep",sim,".asc",sep=""))
    for(rep in 1:10){
      presences <- getSample(1000,calStack,occurrences,"presence")
      background <- getSample(10000,calStack,occurrences,"background")
      paSample <- rbind(presences,background)
      # Generate occupied niche and environmental range models
      #on <- genON(presences,background,calStack,xferStack,pair[1],pair[2],repCount)
      # Calibrate Maxnet SDM
      m1 <- runMaxnetModel(paSample,calStack,xferStack,nicheStack,pair[1],pair[2],repCount)
      m2 <- runRFModel(paSample,calStack,xferStack,nicheStack,pair[1],pair[2],repCount)
      
      # # Evaluate model
      # metricsMAX <- calculateMetrics(m1,on,test,'maxnet',repCount,pair[1],pair[2])
      # metricsRF <- calculateMetrics(m2,on,test,'rf',repCount,pair[1],pair[2])
      # metrics <- rbind(metricsMAX,metricsRF)
      # 
      # if (typeof(results) == "list"){
      #   results <- rbind(results,metrics)
      # }else{
      #   results <- metrics
      #}
      repCount <- repCount + 1
    }
  }
}
write.csv(results,"F:\\PhD\\prediction_extent\\results2.csv")
