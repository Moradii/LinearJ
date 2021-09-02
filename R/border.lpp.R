#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.linnet
#' @import spatstat.core
#' @import stats
#' @import utils
#' @export
border.lpp <- function(X){
  
  if (any(class(X)=="linnet")){
    
    vdx <- vertexdegree(X)
    p <- X$vertices
    l <- X
    
    vers <- as.data.frame(p)
    vers <- vers[which(vdx==1),]
    
    borderpoints <- lpp(vers,l,check=FALSE)
    
  } else if (any(class(X)=="lpp")){
    
    vdx <- vertexdegree(X$domain)
    p <- X$domain$vertices
    l <- domain(X)
    
    vers <- as.data.frame(p)
    vers <- vers[which(vdx==1),]
    
    borderpoints <- lpp(vers,l,check=FALSE)
  }
  else{
    stop("X should be either of class linnet or lpp")
  }
  
  return(borderpoints)
}
