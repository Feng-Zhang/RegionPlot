\name{DurocGff}
\alias{DurocGff}
\docType{data}
\title{
the gene structure information in Duroc
}
\description{
This data set gives the annotated gene structure information in Duroc, which reformulated from Sus scrofa 10.2.73.gtf in Ensemble 
}
\usage{DurocGff}
\format{
  A data frame with 191919 observations on the following 6 variables.
  \describe{
    \item{\code{chr}}{a character vector}
    \item{\code{structure}}{a character vector}
    \item{\code{leftPos}}{a numeric vector}
    \item{\code{rightPos}}{a numeric vector}
    \item{\code{strand}}{a character vector}
    \item{\code{geneName}}{a character vector}
  }
}

\source{\url{http://asia.ensembl.org/index.html}}

\examples{
data(DurocGff)

}
\keyword{datasets}
