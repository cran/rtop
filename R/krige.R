ffun = function(...) {
  rkrige(...)
}


rkrige = function(observations, obs0, obscors, newcor, vObs, c0arr, nmax, inew, cv, 
                  unc0, mdist, maxdist, singMat, varInv, wlim, debug.level,
                  wlimMethod, simul = FALSE, BLUE = FALSE, varClean = FALSE, corlines = NULL, remNeigh = FALSE) {
  
  
  
  naobs = which(is.na(obs0)) 
  naobs = unique(c(naobs, corlines))
  if (varClean & FALSE) {
    vm = vObs
    diag(vm) = 1
    vm[upper.tri(vm)] = 1
    mins = apply(vm, MARGIN = 1, FUN = function(x) min(x))
    naobs = c(naobs, which(mins < 1e-9))
  }
  
  if (length(naobs) > 0) {
    
    if (debug.level > 1) observations = observations[-naobs]
    obs0 = obs0[-naobs]
    obscors = obscors[-naobs,]
    vObs = vObs[-naobs,-naobs]
    c0arr = c0arr[-naobs]
    unc0 = unc0[-naobs]
  }
  nobs = length(obs0)
  nneigh = nobs    
  neigh = c(1:nobs)
  if (!singMat) {
    if (nobs <= nmax && mdist < maxdist) {
      if (cv) {
        vMat = vObs[-inew,-inew]
        c0arr = c0arr[-inew]
        obs = obs0[-inew]
        unc = unc0[-inew]
        neigh = neigh[-inew]
      } else {
        vMat = vObs
        obs = obs0
        unc = unc0
      }
    } else {
      #  There are limits on distance or numbers
      if (mdist > maxdist) {
        distm = spDistsN1(obscors,newcor)
        neigh = which(distm < maxdist)
      }
      if (cv) neigh = neigh[!neigh%in%inew]
      if (nobs > nmax) {
        cOrder = order(c0arr)   
        neigh = cOrder[cOrder %in% neigh][1:nmax]
      }
      neigh = neigh[!is.na(neigh)]
      if (length(neigh) < 1) {
        warning(paste("No neighbours for area", inew, "within maxdist", maxdist))
        return(list(pred = c(NA, krigingError=NA, slambda=NA ),
                    lambda = NA, c0arr = NA, obs = NA, unc = NA, nneigh = NA, 
                    neigh = NA, mu = BLUE))
      }
      
      if (length(neigh) <= nobs) {
        c0arr = c0arr[neigh]
        vMat = vObs[neigh,neigh]
        obs = obs0[neigh]
        unc = unc0[neigh]
      }
    }
    nneigh = length(c0arr)
    if (BLUE) vInv = try(solve(vMat))
    vMat = rbind(vMat,1)
    vMat = cbind(vMat,1)
    
    diag(vMat)[1:nneigh] = -unc
    vMat[nneigh+1,nneigh+1] = 0
    
    while (TRUE) {
      varInv = try(solve(vMat))
      if (is(varInv,"try-error") || (BLUE && is(vInv, "try-error"))) {
        emsg = paste("Error in solve.default(vMat) : \n",
                     "system is computationally singular.\n",
                     #                  "Error most likely occured because two or more areas/lines being (almost) identical \n",
                     "Error most likely occured because two or more areas being (almost) identical \n",
                     "checking prediction location", inew, "\n neighbours",paste(neigh, collapse = " "),"\n",
                     "variance matrix", "\n")  
        for (irr in 1:dim(vMat)[1]) emsg = paste(emsg, paste(vMat[irr,], collapse = " "), "\n")
        stop(emsg)                              
      }
      if (!remNeigh) break
      
      lambda = varInv %*% c(c0arr, 1)
      slambda = sum(abs(lambda[1:(length(lambda)-1)]))
      nneigh = length(c0arr)
      if (slambda < wlim) break
      
      vm = vMat[1:nneigh, 1:nneigh]
      for (ii in 1:nneigh) for (jj in 1:nneigh) vm[ii,jj] = cor(vMat[ii,1:nneigh], vMat[jj, 1:nneigh])
#      vm = cor(vMat)
      diag(vm) = 0
      vm[upper.tri(vm)] = 0
      maxs = apply(vm, MARGIN = 1, FUN = function(x) max(x))
      if (max(maxs) < 0.9) break
      rn = which.max(maxs)
      vMat = vMat[-rn, -rn]
      c0arr = c0arr[-rn]
      obs = obs[-rn]
    }
  } else {
    vMat = vObs
    obs = obs0
    unc = unc0
  }
  c0arr[nneigh+1] = 1
  lambda = varInv %*% c0arr
  slambda = sum(abs(lambda[1:nneigh]))
  if (BLUE) BLUE = sum(vInv %*% c0arr)/sum(vInv)
  s0lambda = slambda
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
  krigingError = sum(lambda*c0arr) 
  
  if (debug.level >1) {
    distm = spDistsN1(obscors,newcor)[neigh]
    if (simul) {
      lobs = obs
    } else {
      lobs = observations[neigh,]
    }
    lobs = rbind(lobs, mu =  rep(0, (dim(lobs)[2])))
    lobs = cbind(lobs, data.frame(id = c(neigh, 0), 
                                  edist = c(distm, 0), lambda = lambda, c0 = c0arr, 
                                  obs = c(obs, 1), unc = c(unc, 0), 
                                  lambda_times_obs = lambda*c(obs, 0))) 
    print("neighbours")
    print(lobs, 3)
    print("covariance matrix ")
    print(vMat,3)     
  }
  
  list(pred = c(sum(lambda[1:nneigh] * obs), krigingError, slambda ),
       lambda = lambda, c0arr = c0arr, obs = obs, unc = unc, nneigh = nneigh, 
       neigh = neigh, mu = BLUE)
}