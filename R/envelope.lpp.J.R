#' @export
J.envelope.lpp <- function(X,Y,lambda=NULL,fun=LinearJinhom,nsim=19,nrank=1,Interval=TRUE,prob=0.025,seed=1234,verbose=TRUE,...){

  L <- as.linnet(X)


  if(missing(Y)){

    if(is.null(lambda)) lambda <- npoints(X)/volume(X$domain)

    set.seed(seed)
    Y <- rpoislpp(lambda = lambda,L=L,nsim = nsim)
    # x <- replicate(nsim, rpoislpp(lambda, L), simplify=FALSE)
  }
  else{
    nsim <- length(Y)
  }


  out <- vector("list", nsim)
  M <- matrix(nrow =513 ,ncol = nsim)

  lower <- integer()
  upper <- integer()

  if(missing(fun)) fun <- LinearJinhom

  JX <- fun(X,lambda=NULL,...)
  r <- attr(JX,"r")

  M <- lapply(X=1:nsim, function(i){

    if(verbose){
      if(i<nsim) {
        cat(paste(i),",")
        flush.console()
      } else {
        cat(paste(i),"\n")
        flush.console()
      }}

    fun(Y[[i]],r=r,lambda=NULL,...)$LinearJinhom
  })

  M <- do.call(cbind,M)

  if(Interval){
    for (j in 1:nrow(M)){
      row <- sort(M[j,])
      lower[j] <- quantile(row,probs = prob)
      upper[j] <- quantile(row,probs = 1-prob)
    }
  }
  else{
    for (j in 1:nrow(M)){
      row <- sort(M[j,])
      lower[j] <- row[nrank]
      upper[j] <- row[nsim-nrank+1]
    }
  }

  out <- list(lower=lower,Upper=upper,r=r)

  if(Interval) attr(out,"prob") <- prob

  attr(out,"Jvalues") <- M
  attr(out,"Jmain") <- JX
  attr(out,"sims") <- Y

  if(!Interval) attr(out,"nrank") <- nrank

  class(out) <- "envelopeJ"
  return(out)
}

print.envelopeJ <- function(X){
  cat("envelope for inohomogeneous Linear J-function");
}
#' @export
plot.envelopeJ <- function(X,cex.legend=1,...){

  plot(X$r,X$lower,type = "n",xlab = "r",
       ylab=expression(italic(J[L][","][inhom](r))),
       ylim = c(min(X$lower,attr(X,"Jmain")$LinearJinhom),max(X$Upper,attr(X,"Jmain")$LinearJinhom)),
       ...)

  lines(X$r,X$lower,type="l",col="white")
  lines(X$r,X$Upper,type="l",col="white")

  polygon(c(X$r, rev(X$r)), c(X$Upper, rev(X$lower)),col = "grey70", border = NA)

  points(X$r,attr(X,"Jmain")$LinearJinhom,type="l",...)
  points(X$r,rep(1,length(X$r)),type="l",lty=2,col=2,...)

  if (abs(max(X$Upper,attr(X,"Jmain")$LinearJinhom)-1) > abs(min(X$lower,attr(X,"Jmain")$LinearJinhom)-1)){

    legend("topleft",inset = 0.05,
           legend = c(expression(italic(hat(J)[L][","][inhom](r))),
                      expression(italic({J[L][","][inhom]^{theo}}(r))),
                      expression(italic({hat(J)[L][","][inhom]^{high}}(r))),
                      expression(italic({hat(J)[L][","][inhom]^{low}}(r)))),
           col = c(1,2,"grey70","grey70"),lty = c(1,2,1,1),
           pt.cex=1,
           cex=cex.legend)
  }
  else{

    legend("bottomleft",inset = 0.05,
           legend = c(expression(italic(hat(J)[L][","][inhom](r))),
                      expression(italic({J[L][","][inhom]^{theo}}(r))),
                      expression(italic({J[L][","][inhom]^{high}}(r))),
                      expression(italic({J[L][","][inhom]^{low}}(r)))
           )
           ,col = c(1,2,"grey70","grey70"),lty = c(1,2,1,1),
           pt.cex=1,
           cex=cex.legend)

   }
}
