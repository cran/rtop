rtopCluster = function(nclus, ..., action = "start", type, outfile = NULL ) {
  cl = getOption("rtopCluster")
  if (length(cl) > 0 && (action == "stop" | action == "restart")) {
    parallel::stopCluster(cl)
    options(rtopCluster = NULL)
  } 
  if (length(cl) > 0 && action == "start") {
    if (length(list(...)) > 0) parallel::clusterEvalQ(cl, ...)
  } else if (action == "start" | action == "restart") {
    if (!requireNamespace("parallel")) stop("Not able to start cluster, parallel not available")
    if (missing(type) || is.null(type)) {
      cl <- parallel::makeCluster(nclus, outfile = outfile) 
    } else {
      cl <- parallel::makeCluster(nclus, type, outfile = outfile)
    }
#    doParallel::registerDoParallel(cl, nclus)
    if (length(list(...)) > 0) parallel::clusterEvalQ(cl, ...)
    options(rtopCluster = cl)    
  }
  getOption("rtopCluster")
}

