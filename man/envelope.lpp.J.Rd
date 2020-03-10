\name{J.envelope.lpp}
\alias{J.envelope.lpp}
%- Also NEED an '\alias' for EACH other topic documented here.
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
interval.
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
  arguments passed to \code{\link{LinearJinhom}.
}
}
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
