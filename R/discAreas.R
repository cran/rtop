

discBinAreas = function(object,object2,dist,resol,stype) {
  ad = sqrt(object)/2
  ad[2] = sqrt(object2)/2
  Srl = list()
  dAreas = list()
  for (i in 1:2) {
    pt1 = c(0,ifelse(i==1,0,dist))
    x1 = pt1[1]-ad[i]
    x2 = pt1[1]+ad[i]
    y1 = pt1[2]-ad[i]
    y2 = pt1[2]+ad[i]
    boun = data.frame(x=c(x1,x2,x2,x1,x1),y=c(y1,y1,y2,y2,y1))
    boun = Polygon(SpatialPoints(boun))
    dAreas[[i]] = spsample(boun,resol,stype,offset = c(0.5,0.5))
  }
  dAreas
}



rtopDisc.rtopVariogram = function(object, params = list(), ...) {
  params = getRtopParams(params, ...)
  resol = params$hresol ^2
  hstype = params$hstype
# Discretize binned areas from the variogram for pdf or Ghosh calculation
#  rta = list()
#  for (i in 1:dim(object)[1]) {
#    rta[[i]] = rtopDiscAreas(object$a1[i],
#            object$a2[i],object$dist[i],
#            resol = resol,stype = stype)
#a1 = object$a1[i]
#a2 = object$a2[i]
#dist = object$dist[i]
#rtopDiscAreas(a1,a2,dist,resol = resol,stype = stype)
#  }

  mapply(discBinAreas,as.list(object$a1),
            as.list(object$a2),as.list(object$dist),
            MoreArgs = list(resol = resol,stype = hstype),SIMPLIFY = FALSE)
}




rtopDisc.rtop = function(object,params = list(), ...) {
  object$params = getRtopParams(object$params,params, ...)
  observations = object$observations
  if ("predictionLocations" %in% names(object)){
    predictionLocations = object$predictionLocations
    bbo = data.frame(t(bbox(observations)))
    bbp = data.frame(t(bbox(predictionLocations)))
    bb = rbind(bbo,bbp)
  } else bb = bbox(observations)
  coordinates(bb) =  as.formula(paste("~",names(bb)[1],"+",names(bb)[2]))
  object$dObs = rtopDisc(observations,bbox(bb),params = object$params)
  object@observations@data$ddim = unlist(lapply(object$dObs,FUN = function(are) dim(coordinates(are)[1])))
  if ("predictionLocations" %in% names(object)){
    object$dPred = rtopDisc(predictionLocations,bbox(bb),params = object$params)
    object@predictionLocations@data$ddim = unlist(lapply(object$dPred,FUN = function(are) dim(coordinates(are)[1])))
  }
  object
}







rtopDisc.SpatialPolygonsDataFrame = function(object, params = list(), bb = bbox(object), ...) {
  rtopDisc(as(object,"SpatialPolygons"), params = params, bb, ...)
}

rtopDisc.SpatialPolygons = function(object, params = list(), bb = bbox(object), ...) {
  params = getRtopParams(params, ...)
  stype = params$rstype
  resol = params$rresol
  if (stype == "random" | stype == "regular") {
    lapply(object@polygons,FUN=function(pol) spsample(pol,resol,stype,offset=c(0.5,0.5)))
  } else if (stype == "rtop") {
    bbdia = sqrt(bbArea(bb))
    small = bbdia/100
    ires0 = 1
    nps = length(object@polygons)
    spp = vector("list",nps)
    for (ip in 1:nps) {
      lpoly = SpatialPolygons(list(object@polygons[[ip]]))
      ba = bbox(lpoly)
      ipts = resol-1
      ires = ires0
      while (ipts < resol) {
        ires = ires*2
        xd = bbdia/(ires)
        if (bbArea(ba)/(xd*xd) > (resol-2)) {
          x = seq(bb[[1]]-small,bb[[3]]+small,xd)
          y = seq(bb[[2]]-small,bb[[4]]+small,xd)
          x = x[x > ba[[1]] & x < ba[[3]] ]
          y = y[y > ba[[2]] & y < ba[[4]] ]
          pts = expand.grid(x=x,y=y)
          if (dim(pts)[1] >= 1) {
            coordinates(pts) = ~x+y
            pts = pts[!is.na(sp:::pointsInSpatialPolygons(pts,lpoly)),]
#            pts = pts[!is.na(pointsInSpatialPolygons_loc(pts,lpoly)),]
            ipts = dim(coordinates(pts))[1]
          }
        }
      }
      print(paste("Sampling from area number",ip,"containing",ipts,"points"))
      spp[[ip]] = pts
    }
    spp
  } else stop(paste("Unknown sampling type:",stype))
}


#pointsInSpatialPolygons_loc <- function(pts, SpPolygons) {
#    pls = slot(SpPolygons, "polygons")
#    lb <- lapply(pls, function(x) as.double(bbox(x)))
#    cpts <- coordinates(pts)
#    storage.mode(cpts) <- "double"
#    mode.checked <- storage.mode(cpts) == "double"
#    cand0 <- .Call("pointsInBox", lb, cpts[,1], cpts[,2], PACKAGE="sp")
#    m <- length(pls)
#    cand <- .Call("tList", cand0, as.integer(m), PACKAGE="sp")
#    res <- sp:::pointsInPolys2(pls, cand, cpts, mode.checked=mode.checked)
#    res
#}
#

#pointsInSpatialPolygons_local <- function(pts, SpPolygons, pls, lb) {
##    pls = slot(SpPolygons, "polygons")
##    lb <- lapply(pls, function(x) as.double(bbox(x)))
#    cpts <- coordinates(pts)
#    storage.mode(cpts) <- "double"
#    mode.checked <- storage.mode(cpts) == "double"
#    cand0 <- .Call("pointsInBox", lb, cpts[,1], cpts[,2], PACKAGE="sp")
#    m <- length(pls)
#    cand <- .Call("tList", cand0, as.integer(m), PACKAGE="sp")
#    res <- sp:::pointsInPolys2(pls, cand, cpts, mode.checked=mode.checked)
#    res
#}
