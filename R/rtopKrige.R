rtopKrige.rtop = function(object, varMatUpdate = FALSE, ...) {
  observations = object$observations
  predictionLocations = object$predictionLocations
  params = getRtopParams(object$params, ...)
  if (!all(c("varMatObs", "varMatPredObs") %in% names(object)) | varMatUpdate) 
    object = varMat(object, varMatUpdate, ...)
  
  varMatObs = object$varMatObs
  varMatPredObs = object$varMatPredObs
  
  krigeRes = rtopKrige(object = observations, predictionLocations = predictionLocations, 
    varMatObs = varMatObs, 
    varMatPredObs = varMatPredObs, params = params, 
    formulaString = object$formulaString, ...)
  object$predictions = krigeRes$predictions
  if ("cvInfo" %in% names(krigeRes)) object$cvInfo = krigeRes$cvInfo
  object
}  


rtopKrige.SpatialPolygonsDataFrame = function(object, predictionLocations = NULL,
    varMatObs, varMatPredObs, varMat, params = list(), formulaString,  
    sel, ...) {
    rtopKrige.default(object, predictionLocations, varMatObs, 
          varMatPredObs, varMat, params, formulaString,  
          sel, ...) 
}



rtopKrige.default = function(object, predictionLocations = NULL,
    varMatObs, varMatPredObs, varMat, params = list(), formulaString,  
    sel, wret = FALSE, ...) {
  params = getRtopParams(params, ...)
#
  cv = params$cv  
#  else object$params$cv = params$cv = cv
  nmax = params$nmax    
  wlim = params$wlim
  wlimMethod = params$wlimMethod
  maxdist = params$maxdist
  debug.level = params$debug.level
  if (!missing(varMat)) {
    if (missing(varMatObs)) varMatObs = varMat$varMatObs
    if (missing(varMatPredObs)) varMatPredObs = varMat$varMatPredObs
  }
  depVar = as.character(formulaString[[2]])
  observations = object
  obs0 = observations[[depVar]]
  nobs = dim(coordinates(object))[1]
  obscors = coordinates(observations)
  if (cv) newcors = obscors else newcors = coordinates(predictionLocations)

  npred = ifelse(cv,nobs,ifelse(!missing(sel),length(sel),dim(newcors)[1]))
  if (missing(sel)) sel = c(1:npred)
  if (params$unc && "unc" %in% names(observations)) {
    unc0 = observations$unc
  } else unc0 = array(0,nobs)
#  
  mdist = sqrt(bbArea(bbox(observations)))
  if (nobs < nmax && mdist < maxdist && !cv) {
    varMat = rbind(varMatObs,1)
    diag(varMat) = unc0
    varMat = cbind(varMat,1)
    varMat[nobs+1,nobs+1] = 0    
    varInv = solve(varMat)    
    singMat = TRUE 
  } else singMat = FALSE
#
  if (cv) {
    predictionLocations = observations
    varMatPredObs = varMatObs
  }
#  
  predictions = data.frame(var1.pred = rep(0,npred),var1.var = 0,sumWeights = 0)
  if (cv) {
    predictions = cbind(predictions,observed = observations[[depVar]], residual=0, zscore = 0)
    cvInfo0 = data.frame(inew = c(0),obs=c(1),pred=c(0),var=c(0),residual = c(0), zscore = c(0),neigh = c(0),neighobs=c(0),weight=c(0),sumWeight=c(0),c0arr=c(0))
    cvInfo = list()
  }
#
  if (wret) weight = matrix(0,nrow = npred,ncol = nobs)
  for (inew in sel) {         
    if (debug.level > 1) print("\n")
#  for (inew in 1:20) {         
    if (cv) {
      if (debug.level >=1) print(paste("Cross-validating location", inew, 
           " out of ",npred," observation locations"))
      if (debug.level > 1) print(observations@data[inew,] )
#      if (cv == inew && inew > 1) browser()
    } else {
      if (debug.level >=1) print(paste("Predicting location ",inew,
           " out of ", npred," prediction locations" ))
      if (debug.level > 1 && is(predictionLocations, "SpatialPolygonsDataFrame")) 
                print(predictionLocations@data[inew,] )
    }
    newcor = newcors[inew,]
    c0arr = varMatPredObs[,inew]
    nneigh = nobs
    obs = obs0
    unc = unc0
    neigh = c(1:nobs)
    if (!singMat) {
      if (nobs < nmax && mdist < maxdist) {
#  cross-validation, but no limits on distance or numbers
        varMat = varMatObs[-inew,-inew]
        c0arr = c0arr[-inew]
        obs = obs0[-inew]
        unc = unc0[-inew]
        neigh = neigh[-inew]
      } else {
#  There are limits on distance or numbers
        if (mdist > maxdist) {
          distm = spDistsN1(obscors,newcor)
          neigh = which(distm < maxdist)
        }
        if (cv) neigh = neigh[-inew]
        if (nobs > nmax) {
          cOrder = order(c0arr)   
          neigh = cOrder[cOrder %in% neigh][1:nmax]
        }
        if (length(neigh) < nobs) {
          c0arr = c0arr[neigh]
          varMat = varMatObs[neigh,neigh]
          obs = obs0[neigh]
          unc = unc0[neigh]
        }
      }
      nneigh = length(c0arr)
      varMat = rbind(varMat,1)
      varMat = cbind(varMat,1)

      diag(varMat)[1:nneigh] = unc
      varMat[nneigh+1,nneigh+1] = 0
      
      varInv = try(solve(varMat))
      if (is(varInv,"try-error")) stop(paste("Error in solve.default(varMat) : \n",
                  "system is computationally singular.\n",
#                  "Error most likely occured because two or more areas/lines being (almost) identical \n",
                  "Error most likely occured because two or more areas being (almost) identical \n",
                  "checking prediction location", inew, "\n neighbours",paste(neigh, collapse = " ")))
    }
    c0arr[nneigh+1] = 1
    lambda = varInv %*% c0arr
    krigingError = sum(lambda*c0arr) 
    slambda = sum(abs(lambda[1:nneigh]))
    while (slambda > wlim) {
      if (wlimMethod == "all") {
        oslambda = slambda
        lambda = lambda/slambda*(wlim/1.01)
#      lambda[1:nneigh] = lambda[1:nneigh]/slambda*(wlim/1.01)
        lambda[1:nneigh] = lambda[1:nneigh]+(1-sum(lambda[1:nneigh]))/nneigh
        slambda = sum(abs(lambda[1:nneigh]))
      } else if (wlimMethod == "neg") {
        wdiv = 1.1
        oslambda = slambda
        neg =  which(lambda[1:nneigh] < 0)
        pos = which(lambda[1:nneigh] > 0)
        lambda[neg] = lambda[neg] /wdiv
        ndiff = (wdiv-1)*sum(abs(lambda[neg]))
        lambda[pos] = lambda[pos] - ndiff*lambda[pos]/sum(lambda[pos])
        slambda = sum(abs(lambda[1:nneigh]))/1.00001             
     }
      if (debug.level >1) 
         print(paste("optimizing lambdas",oslambda, slambda,sum(lambda[1:nneigh]),lambda[nneigh+1]))
    }
    predictions$var1.pred[inew] = sum(lambda[1:nneigh] * obs)
    predictions$var1.var[inew] = krigingError
    predictions$sumWeights[inew] = slambda
    if (wret) weight[inew,neigh] = lambda[1:nneigh]    
    if (debug.level >1) {
      distm = spDistsN1(obscors,newcor)[neigh]
      lobs = observations@data[neigh,]
      lobs = rbind(lobs, mu =  rep(0, (dim(lobs)[2])))
      lobs = cbind(lobs, data.frame(id = c(neigh, 0), 
                   edist = c(distm, 0), lambda = lambda, c0 = c0arr, 
                  obs = c(obs, 1), unc = c(unc, 0), 
                  lambda_times_obs = lambda*c(obs, 0))) 
      print("neighbours")
      print(lobs, 3)
      print("covariance matrix ")
      print(varMat,3)     
    }
    if (cv) {
      predictions$residual[inew] = observations[[depVar]][inew]-predictions$var1.pred[inew]
      predictions$zscore[inew] = predictions$residual[inew]/sqrt(predictions$var1.var[inew])
      cvNew = data.frame(predictions[inew,],c0arr=lambda[nneigh+1],neigh=0)
      for (jnew in 1:nneigh) cvNew = rbind(cvNew,
          data.frame(var1.pred = 0,
            var1.var = unc[neigh[jnew]], sumWeights=lambda[jnew],
            observed = observations@data[[depVar]][neigh[jnew]], residual = 0, zscore = 0,
            c0arr=c0arr[jnew],neigh=neigh[jnew]))
      if (inew == 1) cvInfo = cvNew else cvInfo = rbind(cvInfo,cvNew)
      if (debug.level > 1) {
        print("prediction") 
        print(cbind(predictionLocations@data[inew,],predictions[inew,]))
      }
    }
  }  
  if ("data" %in% names(getSlots(class(predictionLocations)))) {
    predictionLocations@data = cbind(predictionLocations@data, predictions)
    predictions = predictionLocations
  } else {
    predictions = addAttrToGeom(predictionLocations, predictions, match.ID = FALSE)
  }
  if (wret) {
    return(weight)
  } else if (cv) {
    predictions$observed = observations[[depVar]]
    return(list(predictions = predictions, cvInfo = cvInfo))
  } else return(list(predictions = predictions))
}









      