
library(raster)
library(virtualspecies)
output <- "D:\\PHDExperimentOutputs\\Transferability\\landscapes\\"


for(i in 1:888){
  e1 <- raster(paste(output,'env\\temp_ls',i,".tif",sep=""))
  e2 <- raster(paste(output,'env\\precip_ls',i,".tif",sep=""))
  e3 <- raster(paste(output,'env\\habitat_ls',i,".tif",sep=""))
  l1_predictors <- stack(e1,e2,e3)
  names(l1_predictors) <- c('temp','precip','habitat')
  # Formatting of the response functions
  species1.parameters <- formatFunctions(temp = c(fun = 'dnorm', mean = 0.65, sd = 0.3),
                                         precip = c(fun = 'dnorm', mean = 0.35, sd = 0.3),
                                         habitat = c(fun = 'dnorm', mean = 0.65, sd = 0.3))
  # Generation of the virtual species
  species1 <- generateSpFromFun(raster.stack = l1_predictors,parameters = species1.parameters,plot=FALSE,rescale=FALSE)
  writeRaster(species1$suitab.raster,paste(output,'suitability\\suitability',i,'.asc',sep=''),format='ascii',overwrite=TRUE)
}



                              