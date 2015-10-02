\name{colorhcplot}
\alias{colorhcplot}
\title{
Colorful Hierarchical Clustering Dendrograms
}
\description{
This function takes a "hclust-class" object and a factor describing the groups as arguments and generates colorful dendrograms in which leaves belonging to different groups are identified by colors. This function produces a plot that allows to easily detect if leaves that clustered together also belong to the same group.
}
\usage{
colorhcplot(hc, fac, hang = 0.1, main = "Cluster Dendrogram", color = T, 
            lab.cex = 0.9, ax.low = "auto", ax.top = "auto", lwd = 3, 
            window.adj = 1, tick.at = 5)
}

\arguments{
  \item{hc}{
is an object of class hclust, typically the result of a hclust() function call.
}
  \item{fac}{
is a factor and defines the grouping.
}
  \item{hang}{
hang as in \code{\link{hclust}}. hang is the fraction of the plot height by which labels should hang below the rest of the plot. A negative value will cause the labels to hang down from 0.
}
  \item{main}{
title of the dendrogram plot. Default value is "Cluster Dendrogram".
}
  \item{color}{
logical whose default value is TRUE. When color is set to FALSE, only leaf labels will be colored according to the grouping, while leaves/branches will be black.
}
  \item{lab.cex}{
numeric value for adjusting the font size of the leaf labels (and legend text).
}
  \item{ax.low}{
numeric value that sets the height of the lower end of the plot y-axis. Alternatively, when ax.low is set to "auto" (default value), the scale is set automatically.
}
  \item{ax.top}{
numeric value that sets the height of the upper end of the plot y-axis (as for ax.low). 
}
  \item{lwd}{
is a numeric value and defines the width (in points) of the lines of the dendogram.
}
  \item{window.adj}{
is a numeric value and is used to adjust the margins of the dendrogram plot (zoom-in or zoom-out adjustment).
}
  \item{tick.at}{
is a numeric value and sets the spacing of the tick marks on the y-axis of the plot.
}
}
\details{
In order to generate a colorful dendrogram, colorhcplot() function requires 2 mandatory arguments: hc and fac. hc is the result of a hclust() call, while fac is a factor defining the groups. Therefore, the number of leaves of the dendrogram has to be identical to the length of fac (i.e., length(hc$labels) == length(fac) has to be TRUE).
colorhcplot() employs the default R palette. Therefore, only 8 colors are available for the different groups.   
}
\value{
Calling colorhcplot() returns a colorful dendrogram plot
}
\references{
Online colorhcplot() function reference at: 
\url{http://www.biotechworld.it/bioinf/2015/09/30/colorful-hierarchical-clustering-dendrograms-with-r}
}
\author{
Damiano Fantini \email{damiano.fantini@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hclust}}, 
\code{\link{palette}}
}

\examples{
## colorhcplot(hc, fac) -- examples
##
## Example code will generate 2 dendrograms: 
## 1) a simple color dendrogram using default parameters; 
## 2) a y-axis adjusted dendrogram with thin black lines. 
##
data(geneData, package="Biobase")
geneData <- log2(geneData[,7:19])
rmv.id <- apply(geneData,2,(function(x){which(sum(is.na(x))>0)}))
geneData <- geneData[-rmv.id,]
hc <- hclust(dist(t(geneData)))
fac <- factor(c("canc","canc","norm","canc","norm","norm","canc",
                "norm","canc","norm","norm","canc","canc"),
              levels = c("norm","canc"))
par(mfrow=c(1,2))
## first plot
colorhcplot(hc, fac, main ="default")
## second plot
colorhcplot(hc, fac, main ="color=F", hang=0.1, lab.cex = 0.98, 
            ax.low=16, ax.top=36, tick.at=4, 
            window.adj=1.1, lwd =1, color = FALSE)
}

\keyword{hierarchical}
\keyword{cluster}
\keyword{color}
\keyword{dendrogram}
\keyword{plot}
\keyword{multivariate}
