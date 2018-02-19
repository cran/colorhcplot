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
This package contains the function colorhcplot. This function generates simple colorful dendrograms and requires only 2 mandatory arguments: hc and fac. The argument hc is the result of a hclust() call, while fac is a factor defining the groups. Therefore, the number of leaves of the dendrogram has to be identical to the length of fac (i.e., length(hc$labels) == length(fac) has to be TRUE). The function colorhcplot() employs a custom color palette. However, users can specify a custom list of colors.  
}
\author{
\packageAuthor{colorhcplot}

Maintainer: \packageMaintainer{colorhcplot}
}
\seealso{
\code{\link{colorhcplot}}
}
