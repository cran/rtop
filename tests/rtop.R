set.seed(1501)
#-----------------------------
## IGNORE_RDIFF_BEGIN
library(rtop)
if (interactive()) options(error = recover)
  # Read directly from shape-files in data directory
  rpath = system.file("extdata",package="rtop")
  library(sf)
  observations = st_read(rpath, "observations")
  predictionLocations = st_read(rpath, "predictionLocations")
## IGNORE_RDIFF_END
  
  
  observations = as(observations, "Spatial")
  predictionLocations = as(predictionLocations, "Spatial")
#Finding a few prediction locations of them
  
  observations = observations[1:30,]
  predictionLocations = predictionLocations[1:2,]
  
  observations$obs = observations$QSUMMER_OB/observations$AREASQKM
  
  # Setting some parameters 
  params = list(gDist = TRUE, cloud = FALSE, rresol = 25, hresol = 3, debug.level = -1)
  # Build an object
  rtopObj = createRtopObject(observations,predictionLocations, params = params, formulaString = "obs ~ 1" )
  # Fit a variogram (function also creates it)
  rtopObj = rtopFitVariogram(rtopObj, iprint = -1)
  print(rtopObj$variogramModel, 3)
  #rtopObj = checkVario(rtopObj)
  rtopObj2 = rtopKrige(rtopObj, cv = TRUE)


  print(attr(rtopObj2$varMatObs,"variogramModel"), 3)
  
  rtopObj3 = rtopKrige(rtopObj)


 varmat = varMat(observations, predictionLocations, variogramModel = rtopObj$variogramModel, 
                 gDistEst = TRUE, gDistPred = TRUE, rresol = 25, hresol = 3)

all.equal(varmat$varMatObs, rtopObj2$varMatObs)
rtopObj4 = rtopKrige(rtopObj2)

#debug(rtop:::rtopDisc.SpatialPolygons)
#  rtopObj5 = rtopKrige(rtopObj, params = list(cnAreas = 5, cDlim = 10, nclus = 2))
  
  print(summary(rtopObj2$predictions))
  print(summary(rtopObj3$predictions))
  print(summary(rtopObj4$predictions))
  print(all.equal(rtopObj4$predictions, rtopObj3$predictions))
  #spplot(rtopObj$predictions,col.regions = bpy.colors(), c("var1.pred","var1.var"))
  
  # Cross-validation
  #spplot(rtopObj2$predictions,col.regions = bpy.colors(), c("observed","var1.pred"))
  print(cor(rtopObj2$predictions$observed,rtopObj2$predictions$var1.pred))
  



  set.seed(1501)
  library(intamap)
  useRtopWithIntamap()
## IGNORE_RDIFF_BEGIN
  output = interpolate(observations,predictionLocations,
     optList = list(formulaString = obs~1, gDist = TRUE, cloud = FALSE, nmax = 10, rresol = 25, hresol = 3), 
        methodName = "rtop", iprint = -1)
## IGNORE_RDIFF_END
  

  print(all.equal(rtopObj4$predictions@data$var1.pred, output$predictions@data$var1.pred))
  print(all.equal(rtopObj4$predictions@data$var1.var, output$predictions@data$var1.var))


# Updating variogramModel
  
  rtopObj5 = varMat(rtopObj4)
  rtopObj6 = updateRtopVariogram(rtopObj5, exp = 1.5, action = "mult")
  rtopObj7 = varMat(rtopObj6)


  
  
#  observations$obs = log(observations$obs)
  
  # Setting some parameters 
  # Build an object
  rtopObj = createRtopObject(observations,predictionLocations, params = params, formulaString = "obs~1")
  # Fit a variogram (function also creates it)
  rtopObj = rtopFitVariogram(rtopObj, iprint = -1)
  #rtopObj = checkVario(rtopObj)

rtopObj10 = rtopSim(rtopObj, nsim = 5, logdist = TRUE, debug.level = -1)
rtopObj11 = rtopObj
rtopObj11$predictionLocations = rtopObj11$observations
#rtopObj11$observations = NULL
rtopObj11$observations$unc = var(rtopObj10$observations$obs)*min(rtopObj10$observations$area)/rtopObj10$observations$area
rtopObj11$predictionLocations$replaceNumber = 1:dim(rtopObj11$predictionLocations)[1]
rtopObj12 = rtopSim(rtopObj11, nsim = 10, replace = TRUE, debug.level = -1)

print(rtopObj10$simulations@data, digits = 3)
print(rtopObj12$simulations@data, digits = 3)


