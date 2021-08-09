\name{rSSI.lpp}
\alias{rSSI.lpp}
\title{
Simulate Simple Sequential Inhibition on a Linear Network}
\description{
Simulate Simple Sequential Inhibition point pattern on a Linear Network.
}
\usage{
rSSI.lpp(n, r, L, giveup = 1000, nsim = 1, ...)
}
\arguments{
  \item{n}{
number of points.
}
  \item{r}{
inhibition distance.
}
  \item{L}{
linear network.
}
  \item{giveup}{
number of rejected proposals.
}
  \item{nsim}{
number of simulation.
}
\item{...}{
ignored.
}
}

\details{
This function simulates Simple Sequential Inhibition point patterns on a given linear network. The points are at least r units far from each other. Take the relationship of n and r into account with respect to the total length of the underlying linear network.
}

\author{
Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\references{
Cronie, O., Moradi, M., and Mateu, J (2020) Inhomogeneous higher-order summary statistics for point processes on linear networks. Statistics and Computing.
}

\value{
an object of class \link[spatstat.linnet]{lpp}.
}

\seealso{
\link[spatstat.core]{rSSI}
}

\examples{
x <- rSSI.lpp(30,.75,easynet)
plot(x)
}
