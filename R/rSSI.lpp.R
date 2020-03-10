#' @export
rSSI.lpp <-  function(n,r,L,giveup = 1000,nsim=1){

  if(nsim>1){
    return(lapply(X=1:nsim, function(i){
      rSSI.lpp(n=n,r=r,L=L,giveup = giveup,nsim=1)
    }))
  }

  if((n*r)> sum(lengths.psp(L$lines)))
    stop("n should be less than tha ratio of total length of network over r")
  m <- 0
  X <- as.lpp(runifpointOnLines(1,L),L=L)
  while (npoints(X) < n){
    Y <- as.lpp(runifpointOnLines(1,L),L=L)
    d <- crossdist.lpp(X,Y)
    if (min(d[d>0]) >= r){
      X <-  superimpose.lpp(X,Y)
    }
    else{
      m <- m+1
      if(m >= giveup) stop("gave up!")
    }
  }
  return(X)
}
