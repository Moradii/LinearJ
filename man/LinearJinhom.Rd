\name{LinearJinhom}
\alias{LinearJinhom}
\title{
Inhomogeneous Linear J-function for Point Processes on Linear Networks}
\description{
Inhomogeneous Linear J-function for point processes on linear networks.
}
\usage{
LinearJinhom(X, lambda = NULL, densitymethod = c("kernel", "Voronoi"),
              r = NULL, rmax = NULL, f = 0.2, nrep = 200,
              distancetype = c("path", "euclidean"),
              bw = c("scott", "lppl"), ...)
}
\arguments{
  \item{X}{
point pattern on a linear network.
}
  \item{lambda}{
intensity fucnyion.
}
  \item{densitymethod}{
type of method to estimate the intensity.
}
  \item{r}{
distance vector to estimate J-function based on.
}
  \item{rmax}{
max of r.
}
  \item{f}{
retention probability if densitymethod is Voronoi.
}
  \item{nrep}{
number of thinning if densitymethod is Voronoi.}
  \item{distancetype}{
type of distance metric.
}
  \item{bw}{
bandwidth selection method if densitymethod is kernel.}
  \item{\dots}{
argumets passed to density estimation
}
}
\details{
Geometrically corrected inhomogeneous Linear J-function for point processes on linear networks.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
