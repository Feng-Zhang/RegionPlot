\name{WzsGff}
\alias{WzsGff}
\docType{data}
\title{
the gene structure information in WuZhiShan
}
\description{
This data set gives the annotated gene structure information in WuZhiShan
}
\usage{WzsGff}
\format{
  A data frame with 226466 observations on the following 6 variables.
  \describe{
    \item{\code{chr}}{a character vector}
    \item{\code{structure}}{a character vector}
    \item{\code{leftPos}}{a numeric vector}
    \item{\code{rightPos}}{a numeric vector}
    \item{\code{strand}}{a character vector}
    \item{\code{geneName}}{a character vector}
  }
}

\source{
Fang X.D. et al. (2012) GigaScience 1, 16.
}

\examples{
data(WzsGff)
}
\keyword{datasets}
