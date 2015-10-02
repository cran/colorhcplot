\name{colorhcplot-package}
\alias{colorhcplot-package}
\docType{package}
\title{
Introduction to the COLORHCPLOT Package
}
\description{
This is a simple one-function package. Please, refer to the colorhcplot() function manual to check how the function works.
}
\details{
\tabular{ll}{
Package: \tab colorhcplot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-30\cr
License: \tab GPL-2 \cr
}
This package contains the function \code{\link{colorhcplot}}. This function generates simple colorful dendrograms and requires only 2 mandatory arguments: hc and fac. 
The argument hc is the result of a hclust() call, while fac is a factor defining the groups. Therefore, the number of leaves of the dendrogram has to be identical to the length of fac (i.e., length(hc$labels) == length(fac) has to be TRUE).
The function colorhcplot() employs the default R palette. Therefore, only 8 colors are available for the different groups. 
}
\author{

Damiano Fantini \email{damiano.fantini@gmail.com}
}

\keyword{package}
\keyword{colorhcplot}

\seealso{
\code{\link{colorhcplot}}

}


