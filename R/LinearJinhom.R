#' @import spatstat
#' @import stats
#' @import utils
#' @import graphics
#' @export
LinearJinhom <- function(X,lambda=NULL,densitymethod=c("kernel", "Voronoi"),r=NULL,rmax=NULL,f=.2,nrep=200,
                          distancetype=c("path","euclidean"),bw=c("scott","lppl"),...){

  stopifnot(is.lpp(X)) # class X checking

  if(is.null(lambda)){ # intensity estimation

    if(missing(densitymethod)) densitymethod <- "kernel"

    if(densitymethod=="kernel"){
      if(missing(bw)) bw <- "scott"

      if(bw=="scott") Xsigma <- bw.scott.iso(X)
      if(bw=="lppl") Xsigma <- bw.lppl(X)

      lambda <- density.lpp(X,sigma = Xsigma,distance = "euclidean",...)
    }
    else{
      lambda <- densityVoronoi.lpp(X,f=f,nrep = nrep,...)
    }
}



  n <- npoints(X)
  L <- as.linnet(X)     # pixellate linear network


  Llines <- as.psp(L)

  # creating a grid in L
  linemask <- as.mask.psp(Llines,dimyx=64)
  lineimage <- as.im(linemask)

  xx <- raster.x(linemask)
  yy <- raster.y(linemask)

  mm <- linemask$m
  xx <- as.vector(xx[mm])
  yy <- as.vector(yy[mm])

  pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=FALSE)
  pixdf <- data.frame(xc=xx, yc=yy)
  p2s <- project2segment(pixelcentres, Llines)


  gridx <- p2s$Xproj$x
  gridy <- p2s$Xproj$y
  grid  <- lpp(data.frame(gridx,gridy),L)
  grid  <- unique(grid,warn = F)
  ngrid <- npoints(grid)

  # distance measuring
  if(missing(distancetype)) distancetype <- "path"

  if(distancetype=="path"){
    dc <- crossdist.lpp(X,grid)
    dp <- pairdist.lpp(X)
  }

  else{
    dc <- crossdist.ppp(as.ppp(X),as.ppp(grid))
    dp <- pairdist.ppp(as.ppp(X))
  }


  # intensity at data/grid points
  lambdap <- lambda[X]  # intensity at data points
  if(!(length(lambdap)==n) | anyNA(lambdap)) stop("NA intensity value")
  Intgrid <- lambda[grid]  # intensity at grid points

  # \bar{\lambda}
  lambdabar <- min(Intgrid)

  # border extraction
  Xborder <- border.lpp(X)

  # measure distance from data/grid points to border
  if(distancetype=="path"){
    dcrossX <- crossdist.lpp(X,Xborder)
    Xtoborder <- apply(dcrossX, 1, FUN=min)

    dcrossgrid <- crossdist.lpp(grid,Xborder) # distance of grid points to border
    gridtoborder <- apply(dcrossgrid, 1, FUN=min)
  }

  else{
    dcrossX <- crossdist.ppp(as.ppp(X),as.ppp(Xborder))
    Xtoborder <- apply(dcrossX, 1, FUN=min)

    dcrossgrid <- crossdist.ppp(as.ppp(grid),as.ppp(Xborder)) # distance of grid points to border
    gridtoborder <- apply(dcrossgrid, 1, FUN=min)
  }


  # r calculation
  if(is.null(r)) {
    if(is.null(rmax)) rmax <- quantile(Xtoborder,.6)
    r <- seq(0.01,rmax,(rmax-0.01)/512)
  }

  # nearest neighbour calculation
  Hlpp <- parallel::mclapply(X=1:length(r), function(i){
    sum <- 0
    OK <- (Xtoborder > r[i])
    H_seq_OK <- which(OK)

    if (sum(OK)==0) {
      return(0)
    }
    else{
      for (j in H_seq_OK){
        # if (Xtoborder[j] >= r[i]) {
        row <- dp[,j]
        p <- which(0<row & row<=r[i])

        if(length(p)>=1){
          denum <- numeric()
          denum <- unlist(lapply(X=1:length(p), function(k){
            countends(L,X[j],dp[p[k],j])
          }))
          sum <- sum+prod(1-(lambdabar/(lambdap[p]*denum)))
          # print(prod (1-(lambdabar/lambdap[p])))
        }
        else{
          sum <- sum+1
        }
        # }
      }
      return(sum/(sum(OK)))
    }

  })

  Hlpp <- unlist(Hlpp)
  ################################################################
  # empty space calculation
  Flpp <- parallel::mclapply(X=1:length(r), function(i){
    sum <- 0
    OK <- (gridtoborder > r[i])
    F_seq_OK <- which(OK)

    if (sum(OK)==0) {
      return(0)
    }
    else{
      for (j in F_seq_OK){

        # if (gridtoborder[j]>r[i]){

        row <- dc[,j]
        p <- which(0<row & row <= r[i])

        if(length(p)>=1){
          denum <- numeric()
          denum <- unlist(lapply(X=1:length(p), function(k){
            countends(L,grid[j],dc[p[k],j])
          }))
          # for (k in 1:length(p)) {
          #   denum[k] <- countends(L,grid[j],dc[p[k],j])
          # }
          sum <- sum+prod(1-(lambdabar/(lambdap[p]*denum)))
        }
        else{
          # denum <- numeric()
          # sum <- sum+prod(1-(lambdabar/(lambdap[p]*denum)))
          sum <- sum+1
        }


        # }
      }
      return(sum/sum(OK))
    }

  })
  Flpp <- unlist(Flpp)

  # puting num and den together to get J

  LinearJinhom <- Hlpp/Flpp

  out <- list(LinearJinhom=LinearJinhom)

  attr(out,"r") <- r
  attr(out,"LinearFinhom") <- 1-Flpp
  attr(out,"LinearGinhom") <- 1-Hlpp

  class(out) <- c("linearJ")

  return(out)

}

print.linearJ <- function(X){
  cat("Inhomogeneous Linear J-function")
}

#' @export
plot.linearJ <- function(x,...){
  plot(attr(x,"r"),x$LinearJinhom,type = "l",xlab = "r",ylab = expression(italic(hat(J)[L][","][inhom](r))),...)
  points(attr(x,"r"),rep(1,length(attr(x,"r"))),lty=2,col=2,type = "l")
  if (abs(max(x$LinearJinhom)-1) > abs(min(x$LinearJinhom)-1)){
    legend("topleft",inset = 0.05,legend = c(expression(italic(hat(J)[L][","][inhom](r))),expression(italic({J[L][","][inhom]^{theo}}(r))))
           ,col = c(1,2),lty = c(1,2))
  }
  else{
    legend("bottomleft",inset = 0.05,legend = c(expression(italic(hat(J)[L][","][inhom](r))),expression(italic({J[L][","][inhom]^{theo}}(r))))
           ,col = c(1,2),lty = c(1,2))
  }
}
