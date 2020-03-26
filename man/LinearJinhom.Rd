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
the estimated intensity to be given when calculating the inhomogeneous Linear J-function. Input should be of class \link[spatstat]{linim}.
}
  \item{densitymethod}{
if no \code{lambda} is given, the intensity will be internally estimated based on either fast kernel intensity estimator by Rakshit et al. (2019) or resample-smoothed Voronoi intensity estimator by Moradi  et al. (2019). If \code{Voronoi} is selected, then arguments \code{f} and \code{nrep} may be selected accordingly.
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
number of thinning if densitymethod is Voronoi.
}
\item{distancetype}{
the type of distance to be given. Difault is the shortest-path distance.
}
\item{bw}{
bandwidth selection method if densitymethod is kernel. Difault is \link[spatstat]{bw.scott.iso}.
}
\item{\dots}{
argumets passed to density estimation.
}
}

\details{
This function computes the geometrically corrected inhomogeneous Linear J-function for point processes on linear networks.

If no estimated intensity is given, this function internally estimates the intensity of the underlying point pattern based on either the function \link[spatstat]{densityQuick.lpp} (\code{densitymethod = "kernel"}) or the function \link[spatstat]{densityVoronoi.lpp} (\code{densitymethod = "Voronoi"}). If \code{densitymethod = "kernel"} is selected, then bandwidth method can be either selected based on \link[spatstat]{bw.scott.iso} (\code{bw=scott}) or \link[spatstat]{bw.lppl} (\code{bw=lppl}).
}
\author{
Mehdi Moradi \email{m2.moradi@yahoo.com}
}
\value{
The estimated inhomogeneous linear J-function.
}

\note{
The function is currently slow.
}

\references{
Cronie, O., Moradi, M., and Mateu, J (2020) Inhomogeneous higher-order summary statistics for point processes on linear networks. Statistics and Computing.
}

\seealso{
\link[spatstat]{bw.scott.iso},\link[spatstat]{bw.lppl},
\link[spatstat]{densityVoronoi.lpp}, \link[spatstat]{densityQuick.lpp}
}
