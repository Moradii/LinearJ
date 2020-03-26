\name{border.lpp}
\alias{border.lpp}
\title{
Border of a Linear Network.
}
\description{
Border of a Linear Network.
}
\usage{
border.lpp(X)
}
\arguments{
  \item{X}{
  either a linear network or a point pattern on a linear network.
}
}
\details{
For a given point pattern or linear network, this function extracts the border of the linear network in consideration which consists of all nodes with degree one.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\value{
an object of class \link[spatstat]{lpp} for which points represent the border points ofthe underlying network.
}

\references{
Cronie, O., Moradi, M., and Mateu, J (2020) Inhomogeneous higher-order summary statistics for point processes on linear networks. Statistics and Computing.
}

\examples{
x <- rpoislpp(10,easynet)
border.lpp(x)
}

\seealso{
\link[spatstat]{lpp}, \link[spatstat]{rSSI}
}
