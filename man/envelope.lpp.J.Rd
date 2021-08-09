\name{J.envelope.lpp}
\alias{J.envelope.lpp}
\title{
Simulation Envelopes of the Inhomogeneous Linear J-function for Point Processes on Linear Networks
}
\description{
Simulation Envelopes of the Inhomogeneous Linear J-function for Point Processes on Linear Networks
}
\usage{
J.envelope.lpp(X, Y, lambda = NULL, fun = LinearJinhom,
                nsim = 19, nrank = 1, Interval = TRUE,
                prob = 0.025, seed = 1234, verbose = TRUE,
                ...)
}
\arguments{
  \item{X}{
a point pattern on a linear network.
}
  \item{Y}{
a list of point pattern (optional).
}
  \item{lambda}{
intensity function.
}
  \item{fun}{
function to build the envelope for.
}
  \item{nsim}{
number of simulation.
}
  \item{nrank}{
rank.
}
  \item{Interval}{
if TRUE, envelopes will be obtained based on \link[stats]{quantile}.
}
  \item{prob}{
probability to create envelope based on quantiles.
}
  \item{seed}{
  seed.
}
  \item{verbose}{
  if TRUE, it shows the progress.
}
  \item{\dots}{
  arguments passed to \link[LinearJ]{LinearJinhom}.
}
}

\author{
Mehdi Moradi \email{m2.moradi@yahoo.com}
}
\value{
The estimated envelope for the inhomogeneous linear J-function.
}

\details{
For a given point pattern \code{X}, this function calculates an envelope based on \code{nsim} point patterns from a Poisson model. \code{nsim} realisations can be either given by user through argument \code{Y}, or be simulated based on a given intensity function \code{lambda}. The linear J-function will be estimated for all \code{nsim} realisations, and then based on all estimated J-functions, lower and upper bound of the envelope will be obtained. Envelope can be obtained using \code{nrank} in which \code{nrank=1} means using min and max. However, envelope can be also obtained based on \link[stats]{quantile} when \code{Interval=TRUE}.
}

\note{
The function is currently slow.
}

\seealso{
\link[LinearJ]{LinearJinhom}, \link[spatstat.linnet]{envelope.lpp}
}

\references{
Cronie, O., Moradi, M., and Mateu, J (2020) Inhomogeneous higher-order summary statistics for point processes on linear networks. Statistics and Computing.
}
