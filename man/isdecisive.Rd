% File src/library/base/man/isdecisive.Rd
\name{isdecisive}
\alias{isdecisive}
\title{isdecisive}
\description{
  Determines whether the dataset is decisive or not and if not it tries to fix the dataset.
  \code{save}.
}
\usage{
  isdecisive(filename, unrooted=T, fflag=F)
}
\arguments{
  \item{filename}{a csv-formatted file that stores the data matrix (dataset).}
  \item{unrooted}{which type of tree to check for.}
  \item{fflag}{a flag to tell the program to fix the dataset.}
}
\value{A \code{\linkS4class{list}} with two elements, \code{Fixed data matrix} and \code{A string as an answer whether the given dataset is decisive or not.} }
\examples{
  library(decisivatoR)
  ans <- isdecisive(filename=system.file("extdata","SmallTest_Unrooted_NotDecisive.csv",package="decisivatoR"),unrooted=T)
}
\keyword{isdecisive}