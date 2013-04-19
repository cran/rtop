rtopVariogram.SpatialPolygonsDataFrame = function(object, ... ) {
  if (missing(object))  stop("rtopVariogram: Observations are missing")
  obs = object@data
  coordinates(obs) = coordinates(object)
  if ("Shape_Area" %in% names(object)) {
    obs$area = object$Shape_Area
  } else obs$area = unlist(lapply(object@polygons,FUN = function(poly) poly@area))
  rtopVariogram(obs, ...)
}





rtopVariogram.SpatialPointsDataFrame = function(object, formulaString, params=list(), cloud, abins, dbins, ...) {
# If params is intamapParams, they will here be included in rtopParams
if (!inherits(params, "rtopParams"))  params = getRtopParams(params, ...)
# amul refers to the number of areal bins per order of magnitude
# dmul refers to the number of distance bins per order of magnitude
amul = params$amul
dmul = params$dmul
if (missing(cloud)) cloud = params$cloud

observations = object
if (missing(observations)) stop("rtopVariogram: Observations are missing")
if (!("area") %in% names(observations) && !("length") %in% names(observations)) 
  stop("rtopVariogram: Observations do not include area (polygons)")
#  stop("rtopVariogram: Observations do not include area (polygons) or length (lines)")
if (missing(formulaString)) {
  if ("obs" %in% names(observations)) { 
    formulaString = "obs ~ 1" 
  } else if ("value" %in% names(observations)) {
    formulaString = "value ~ 1" 
  } else if (length(names(observations@data)) == 1) {
    formulaString = paste(names(observations@data),"~ 1")      
  } else stop("formulaString is missing and cannot be found from data")
  warning(paste("formulaString missing, using",formulaString))      
}
if (!inherits(formulaString,"formula")) formulaString = as.formula(formulaString)

clvar = variogram(formulaString, observations, cloud = TRUE, ...)
.BigInt = attr(clvar, ".BigInt")
clvar$ord = clvar$np
clvar = as.data.frame(clvar)

clvar$a1 = observations@data$area[clvar$left]
clvar$a2 = observations@data$area[clvar$right]


if (cloud) {
  clvar$acl1 = clvar$left
  clvar$acl2 = clvar$right
  clvar = clvar[,-which(names(clvar) %in% c("left","right"))]
  var3d = clvar
  class(var3d) = c("rtopVariogramCloud","data.frame")
  attr(var3d, ".BigInt") = .BigInt
  var3d$np = 1
} else {
  abins = adfunc(NULL, observations, amul)
  dbins = dfunc(NULL, observations, dmul)
  observations$acl = findInterval(observations$area, abins)
  clvar$acl1 = observations$acl[clvar$left]
  clvar$acl2 = observations$acl[clvar$right]

  ich = which(clvar$acl1 > clvar$acl2)
  acl1c = clvar$acl1
  clvar$acl1[ich] = clvar$acl2[ich]
  clvar$acl2[ich] = acl1c[ich]
  
  clvar$dbin = findInterval(clvar$dist, dbins)
  clvar$np = 1
  varnp = aggregate(list(np = clvar$np),list(acl1 = clvar$acl1,acl2 = clvar$acl2,dbin = clvar$dbin),sum)
  var3d = data.frame(np=varnp$np,aggregate(list(dist=clvar$dist,gamma=clvar$gamma,a1 = clvar$a1, a2=clvar$a2),
                       list(acl1 = clvar$acl1,acl2 = clvar$acl2,dbin = clvar$dbin),mean))
  class(var3d) = c("rtopVariogram","data.frame")
}
var3d
}



# Alternative binning:
#  x <- matrix(rnorm(30000), ncol=3)
#  breaks <- seq(-1, 1, length=5)
#  xints <- data.frame(
#  x1=cut(x[, 1], breaks=breaks),
#  x2=cut(x[, 2], breaks=breaks),
#  x3=cut(x[, 3], breaks=breaks))
#  table(complete.cases(xints))
#  xtabs(~ ., xints)




###############################
rtopVariogram.rtop = function(object,... ) {
  observations = object$observations
  formulaString = object$formulaString
  params = object$params

#calling rtopVariogram.SpatialPolygonsDataFrame
  var3d = rtopVariogram(observations,formulaString,params,...)
  if(inherits(var3d,"rtopVariogramCloud"))object$variogramCloud = var3d else object$variogram = var3d
  object
}


################################
