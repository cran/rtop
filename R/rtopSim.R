rtopSim.rtop = function(object, varMatUpdate = FALSE, beta = NA, largeFirst = TRUE, replace = FALSE, 
                        params = list(), dump = NULL, debug.level, ...) {
  params = getRtopParams(object$params, newPar = params,  ...)
  nmax = params$nmax
  cv = params$cv
  maxdist = params$maxdist
  wlim = params$wlim
  wlimMethod = params$wlimMethod
  dots = list(...)
  if (missing(debug.level)) debug.level = params$debug.level
  varClean = params$varClean
  variogramModel = object$variogramModel
  if (is.null(variogramModel)) stop("Cannot do simulations without a variogram model")
  if (length(object$observations) > 0 &(is.null(object$varMatObs) | varMatUpdate))
    object = varMat(object, varMatUpdate, debug.level = debug.level, ...)
  if (is.null(object$varMatPred) || diff(dim(object$varMatPred)) != 0 | varMatUpdate) {
    if (is.null(object$dPred) & !(params$gDistPred & !is.null(object$gDistPred))) {
      object$dPred = rtopDisc(object$predictionLocations, params = params)
    }
    if (params$gDistPred & is.null(object$gDistPred)) {
      object$gDistPred = gDist(object$dPred, params = params)
      varMatPred = object$varMatPred = varMat(object$gDistPred, params = params, 
                                              variogramModel = variogramModel, debug.level = debug.level)
    } else {
      varMatPred = object$varMatPred = varMat(object$dPred, params = params, 
                                              variogramModel = variogramModel, debug.level = debug.level)
    }
  }
  varMatPredObs = object$varMatPredObs
  varMatObs = object$varMatObs
  varMatPred = object$varMatPred
  predictions = object$predictionLocations
  observations = object$observations
  nobs = length(observations)
  predictionLocations = object$predictionLocations
  predictions = predictionLocations
  if (inherits(predictions, "Spatial")) {
  if (!is(predictions, "SpatialPolygonsDataFrame")) {
    aPred = sapply(slot(predictions, "polygons"), function(i) slot(i, "area"))
    predictions = SpatialPolygonsDataFrame(predictions, data = data.frame(area = aPred ))
  } else if (!"area" %in% names(predictions)) {
    predictions$area = sapply(slot(predictions, "polygons"), function(i) slot(i, "area"))
  }
  } else if (inherits(predictions, "sf")) {
    if (!("area" %in% names(predictions))) predictions$area = set_units(st_area(predictions), NULL)
  } else if (inherits(predictions, "sfc_POLYGON")) {
    predictions = st_sf(predictions, area = set_units(st_area(predictions)), NULL)
  } 

  if (replace & !("replaceNumber" %in% names(predictionLocations))) {
    stop("Cannot replace observations if predictionLocations does not have column with replaceNumber")
  } else if (replace && !(dim(observations)[1] >= max(predictionLocations$replaceNumber, na.rm = TRUE) &
                         min(predictionLocations$replaceNumber, na.rm = TRUE) >= 1)) {
    stop("predictionLocations$replaceNumber does not correspond with the number of observations")
  }
  if (!is.null(dump)) save(object, file = dump)    
  singMat = FALSE
  varInv = NULL
  for (isim in 1:params$nsim) {
    predictions$sim = NA
    if (length(dim(observations)) > 0 && dim(observations)[1] > 0) {
      obsall = data.frame(observations)
      obs = obsall[,as.character(object$formulaString[[2]])]
      if (inherits(observations, "Spatial")) {
      obscors = coordinates(observations)
      } else {
        obscors = suppressWarnings(st_coordinates(st_centroid(observations)))
      }
      #      if (params$unc && "unc" %in% names(observations)) {
      #        unc0 = observations$unc
      #      } else unc0 = array(0,nobs)
      nobs0 = dim(observations)[1]
    } else {
      obsall = NULL
      obs = NULL
      obscors = NULL
      nobs0 = 0
    }
    vPred = varMatPred
    vObs = varMatObs
    corlines = NULL
    if (varClean) {
      vm = vObs
      diag(vm) = 1
      vm[upper.tri(vm)] = 1
      mins = apply(vm, MARGIN = 1, FUN = function(x) min(x))
      corlines  = which(mins < 1e-9)
    }
    
    
    vPredObs = varMatPredObs
    dPred = dim(predictions)[1]
    ips = sample(dPred, dPred)
    if (largeFirst) {
      il = which(ips == order(predictions$area, decreasing = TRUE)[1])
      tmp = ips[il] 
      ips[il] = ips[1]
      ips[1] = tmp
    }
    vpo = 1:dPred
    if (interactive() & debug.level) {
      pb = txtProgressBar(1, length(ips), style = 3)
    }
    print(paste0(isim, ". simulation of ", length(ips), " areas"))
    oind = NULL
    if (params$unc && "unc" %in% names(observations)) {
      unc0 = observations$unc
    } else unc0 = array(0,nobs)
    
    for (ip in 1:length(ips)) {
      inew = ips[ip]
      in2 = which(vpo == inew)
      nobs = length(obs)
      if (interactive() & debug.level) setTxtProgressBar(pb, ip)
      if (inherits(predictionLocations, "Spatial")) {
      newcor = coordinates(predictionLocations[inew,])
      } else {
        newcor = suppressWarnings(st_coordinates(st_centroid(predictionLocations[inew,])))
      }
      if (nobs == 0) {
        if (is.na(beta)) stop("No observations found, beta (expected mean) has to be given")
        if (inherits(predictionLocations, "Spatial")) {
          c0 = varioEx(sqrt(bbArea(bbox(predictionLocations[in2,]))), variogramModel)
        } else {
          c0 = varioEx(sqrt(bbArea(st_bbox(predictionLocations[in2,]))), variogramModel)
        }
        inewvar = varMatPred[inew,inew]
        obs = rnorm(1, beta, c0-inewvar)
        vObs = matrix(inewvar, nrow = 1, ncol = 1)
        vPredObs = vPred[inew, -inew, drop = FALSE]
        unc0 = 0
      } else {
        
        mdist = sqrt(diff(range(obscors[,1]))^2 + diff(range(obscors[,2]))^2)
        wlim0 = wlim
        while (TRUE) {
          wlim0 = wlim0/1.05
#          if (ip > 450) browser()
          ret <- try(rkrige(obsall, obs, obscors, newcor, vObs, vPredObs[,in2, drop = FALSE], nmax, inew, cv, 
                     unc0, mdist, maxdist, singMat, varInv, singularSolve = FALSE, wlim0, debug.level, 
                     wlimMethod, simul = TRUE, varClean = FALSE, corlines = corlines, remNeigh = TRUE), silent = TRUE)
          if (is(ret, "try-error")) print(paste("error in simulation of area number", ip))
          if (wlim0 < 1.05 || (!is(ret, "try-error") && ret$pred[2] > 0)) break
        }
        if (!is(ret, "try-error")) {
          nneigh = ret$nneigh
          lambda = ret$lambda
          neigh = ret$neigh
        
          pred = ret$pred
#          if (logdist) {
#            newval = rlnorm(1, log(pred[1]), sqrt(log(pred[2]/pred[1] + 1))) 
#          } else 
         newval = rnorm(1, pred[1], sqrt(pred[2]))
          if (replace && !is.na(predictionLocations$replaceNumber)[inew] && predictionLocations$replaceNumber[inew] != 0) {
            obs[predictionLocations$replaceNumber[inew]] = newval
            oind = c(oind, predictionLocations$replaceNumber[inew])
          } else {
            obs = c(obs, newval)
            oind = c(oind, length(obs))
          }
#          print(paste(ip, length(ip), length(obs), predictionLocations$replaceNumber[inew], newval))
        } else {
          if (replace && !is.na(predictionLocations$replaceNumber)[inew] && predictionLocations$replaceNumber[inew] != 0) {
            obs[predictionLocations$replaceNumber[inew]] = NA
            oind = c(oind, predictionLocations$replaceNumber[inew])
          } else {
            obs = c(obs, NA)
            oind = c(oind, length(obs))
          }
        }
        if (replace && !is.na(predictionLocations$replaceNumber)[inew] && predictionLocations$replaceNumber[inew] != 0) {
          predictionLocations$replaceNumber[inew]
          unc0[predictionLocations$replaceNumber[inew]] = 0  
          vPredObs = vPredObs[,-in2, drop = FALSE]
        } else {
          unc0 = c(unc0, 0)  
          vObs = rbind(vObs, vPredObs[,in2])
          vObs = cbind(vObs, c(vPredObs[,in2], 0))
          nd = dim(vObs)[1]
          if (varClean & any(vObs[nd, 1:(nd-1)] < 1e-9)) corlines = c(corlines, nd)
          vPredObs = vPredObs[,-in2, drop = FALSE]
          vPredObs = rbind(vPredObs, vPred[in2, -in2])
        }
      }  
      obscors = rbind(obscors, newcor)        
      vPred = vPred[-in2, -in2, drop = FALSE]
      vpo = vpo[-in2]      
    }
    if (interactive() & debug.level) close(pb)
    if (replace) {
      predictions$sim[ips] = obs[oind]
    } else predictions$sim[ips] = obs[(nobs0+1) : length(obs)]
    names(predictions)[dim(predictions)[2]] = paste0("sim", isim)
    if (!is.null(dump)) save(object, file = dump)    
  }
  object$simulations = predictions
  object
}
 

rtopSim.default = function(object = NULL, predictionLocations, 
                           varMatObs, varMatPredObs, varMatPred, 
                           variogramModel, ...) {
object = createRtopObject(object, predictionLocations, ...)
if (!missing(varMatObs)) object$varMatObs = varMatObs
if (!missing(varMatPredObs)) object$varMatPredObs = varMatPredObs
if (!missing(varMatPred)) object$varMatPred = varMatPred
object$variogramModel = variogramModel
rtopSim(object, ...)$simulations
}
