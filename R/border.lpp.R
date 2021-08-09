#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.linnet
#' @import spatstat.core
#' @import stats
#' @import utils
#' @export
border.lpp <- function(X){

  if (any(class(X)=="linnet")){
    m <- X$m
    p <- X$vertices
    l <- X
  }

  if (any(class(X)=="lpp")){
    m <- X$domain$m
    p <- X$domain$vertices
    l <- domain(X)
  }

  rowsum <- numeric()

  for (i in 1:length(m[,1])){
    rowsum[i]=sum(m[i,])
  }

  verticesx <- p$x
  verticesy <- p$y
  borderx <- verticesx[which(rowsum==1)]
  bordery <- verticesy[which(rowsum==1)]

   borderpoints <- lpp(data.frame(borderx,bordery),l,check=FALSE)

  return(borderpoints)
}
