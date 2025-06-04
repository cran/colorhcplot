
# ~~~ colorhcplot v.1.4.1 ~~~
# ~~~ by  Damiano Fantini ~~~




#' Colorful Hierarchical Clustering Dendrograms
#' 
#' Build colorful dendrograms based on a "hclust-class" object and a factor 
#' describing the sample groups. 
#' Leaves belonging to different groups are identified by colors, and the 
#' resulting plot enables detection of pure clusters where all leaves 
#' belong to the same group.
#' 
#' 
#' @param hc hclust-class object, typically the 
#' result of a 'hclust()' function call.
#' @param fac factor that defines the grouping.
#' @param hang hang value, as in \code{hclust}. This is the fraction 
#' of the plot height by which labels should hang below the rest of the plot. 
#' A negative value will align all labels at the bottom of the plot.
#' @param main string, title of the dendrogram plot.
#' @param colors NULL or a character vector of length 1 or having the 
#' same length as the number of levels in fac. 
#' This argument defines the palette for the plot.
#' @param lab.cex numeric value for adjusting the font 
#' size of the leaf labels (and legend text).
#' @param ylim numeric, defines the minimum and maximum 
#' value of the y-axis of the plot.
#' @param lwd numeric value that defines the width (in points) 
#' of the lines of the dendogram.
#' @param las numeric value, graphic parameter for the orientation of 
#' the y-axis tick labels.
#' @param lab.mar numeric value, fraction of the plot area that is 
#' reserved for the labels (at the bottom of the plot).
#' 
#' @details In order to generate a colorful dendrogram, 
#' the colorhcplot() function requires 2 mandatory arguments: hc and fac. 
#' hc is the result of a hclust() call, while fac is a factor defining 
#' the groups. The number of leaves of the dendrogram has 
#' to be identical to the length of fac.
#' 
#' @return Calling colorhcplot() returns a colorful dendrogram plot.
#' 
#' @author Damiano Fantini <damiano.fantini@@gmail.com>
#' 
#' @note Online colorhcplot() function reference at: 
#' \url{http://www.biotechworld.it/bioinf/2015/09/30/colorful-hierarchical-clustering-dendrograms-with-r/}
#'
#' @seealso \link{hclust}
#' 
#' @importFrom graphics par axis mtext segments text legend
#' @importFrom methods is
#'
#' 
#' @examples  
#' ### Example 1, using the USArrests dataset
#' data(USArrests)
#' hc <- hclust(dist(USArrests), "ave")
#' fac <- as.factor(c(rep("group 1", 10), 
#'                    rep("group 2", 10),
#'                    rep("unknown", 30)))
#' plot(hc)
#' colorhcplot(hc, fac)
#' colorhcplot(hc, fac, hang = -1, lab.cex = 0.8)
#' 
#' ### Example 2: use the "ward.D2" algorithm and
#' ### the UScitiesD dataset
#' data(UScitiesD)
#' hcity.D2 <- hclust(UScitiesD, "ward.D2")
#' fac.D2 <-as.factor(c(rep("group1", 3), 
#'                      rep("group2", 7)))
#' plot(hcity.D2, hang=-1)
#' colorhcplot(hcity.D2, fac.D2, color = c("chartreuse2", "orange2"))
#' colorhcplot(hcity.D2, fac.D2, color = "gray30", 
#'             lab.cex = 1.2, lab.mar = 0.75)
#'  
#' ### Example 3: use gene expression data 
#' data(geneData, package="colorhcplot")
#' exprs <- geneData$exprs
#' fac <- geneData$fac
#' hc <- hclust(dist(t(exprs)))
#' colorhcplot(hc, fac, main ="default", col = "gray10")
#' colorhcplot(hc, fac, main="Control vs. Tumor Samples") 
#'  
#'  
#'  
#' @export
colorhcplot <- function (hc, fac, hang = 0.1, 
                         main = "Cluster Dendrogram", 
                         colors = NULL, 
                         lab.cex = 1, ylim = NULL, 
                         lwd = 3, las = 1,
                         lab.mar = 0.55) {
  
  stopifnot(methods::is(hc, 'hclust'),
            is.factor(fac))

  # Prepare colors
  all.colPalette <- c("#1f78b4", "#33a02c", "#e31a1c",
                      "#ff7f00", "#6a3d9a", "#b15928",
                      "#a6cee3", "#b2df8a", "#fb9a99",
                      "#fdbf6f", "#cab2d6", "#ffff99")
  
  colLen <- length(levels(fac))
  myNum <- 1 + ceiling((colLen / length(all.colPalette)))
  all.colPalette <- rep(all.colPalette, myNum)       
  
  if(is.null(colors)) {
    fullColors <- all.colPalette[1:length(levels(fac))]
    fullColors <- c("gray70", fullColors)
  } else if (length(colors) == 1){
    fullColors <- all.colPalette[1:length(levels(fac))]
    fullColors <- c("gray70", fullColors)
  } else if(length(colors) ==  length(levels(fac))){
    # good to go
    fullColors <- c("gray70", colors)
  } else {
    stop("The number of colors do not match with the number of levels")
  }
  
  # Reconstruct the dendrogram matrix
  m_ord <- -1 * hc$order
  mgs <- as.vector(hc$merge)
  for (i in 1:length(mgs)) {
    if (mgs[i] < 0) {
      mgs[i] <- as.numeric(-1 * which(m_ord == mgs[i]))
    }
  }
  mgs <- matrix(mgs, ncol = 2)
  mgs <- cbind(mgs, rep(NA, nrow(mgs)))
  for (step in 1:nrow(mgs)) {
    if (mgs[step, 1] < 0 & mgs[step, 2] < 0) {
      mgs[step, 3] <- (-1) * mean(mgs[step, 1:2])
    } else if (mgs[step, 1] * mgs[step, 2] < 0) {
      min <- which(mgs[step, 1:2] < 0)
      plu <- which(mgs[step, 1:2] > 0)
      mgs[step, 3] <- ((-1 * mgs[step, min]) + mgs[mgs[step, 
                                                       plu], 3])/2
    } else {
      mgs[step, 3] <- (mgs[mgs[step, 1], 3] + mgs[mgs[step, 
                                                      2], 3])/2
    }
  }
  mgs <- mgs[, 1:3]
  mgs <- cbind(mgs, matrix(NA, nrow = nrow(mgs), ncol = 3))
  for (step in 1:nrow(mgs)) {
    for (i in 1:2) {
      if (mgs[step, i] < 0) {
        mgs[step, (i + 3)] <- as.numeric(fac[hc$order][(-1) * 
                                                         mgs[step, i]])
      } else {
        mgs[step, (i + 3)] <- mgs[mgs[step, i], 6]
      }
    }
    if (mgs[step, 4] == mgs[step, 5]) {
      mgs[step, 6] <- mgs[step, 5]
    } else {
      mgs[step, 6] <- 0
    }
  }
  mgs[, 6] <- 1 + mgs[, 6]
  dndr_gram <- matrix(NA, ncol = 3, nrow = nrow(mgs))
  colnames(dndr_gram) <- c("x0", "x1", "height")
  dndr_gram[, 3] <- hc$height
  for (step in 1:nrow(mgs)) {
    if (mgs[step, 1] < 0 & mgs[step, 2] < 0) {
      dndr_gram[step, 1] <- mgs[step, 1] * (-1)
      dndr_gram[step, 2] <- mgs[step, 2] * (-1)
    } else if (mgs[step, 1] * mgs[step, 2] < 0) {
      min <- which(mgs[step, 1:2] < 0)
      dndr_gram[step, min] <- mgs[step, min] * (-1)
      plu <- which(mgs[step, 1:2] > 0)
      dndr_gram[step, plu] <- mgs[mgs[step, plu], 3]
    } else {
      dndr_gram[step, 1] <- mgs[mgs[step, 1], 3]
      dndr_gram[step, 2] <- mgs[mgs[step, 2], 3]
    }
  }
  if (hang <= 0) {
    # make all lines end at x-axis
    lv_len <- 0.05 * (max(hc$height) - min(hc$height))
    
  } else {
    lv_len <- hang * (max(hc$height) - min(hc$height))
  }
  
  # Get current settings
  old.mar <- graphics::par()$mar
  graphics::par(mar = c(2.1, 4.1, 4.1, 2.1))
  
  # Define auto.ylims
  auto.ylim <- c(min(dndr_gram[,3] - lv_len), max(dndr_gram[,3]))    
  
  full.range <- auto.ylim[2] - auto.ylim[1]
  if (is.null(ylim)) {
    new.auto.ylim <- c(auto.ylim[1] - (full.range * lab.mar), auto.ylim[2])    
  } else if (is.numeric(ylim) & length(ylim) == 2) {
    new.auto.ylim <- ylim
  } else {
    new.auto.ylim <- c(auto.ylim[1] - (full.range * lab.mar), auto.ylim[2])
  }
  
  # Start plotting
  plot(0:(nrow(mgs) + 2), 0:(nrow(mgs) + 2), xlab = "", 
       ylab = "", 
       type = "n", main = main, ylim = new.auto.ylim, axes = F)
  
  ax2 <- graphics::axis(2, pos = -10000)
  graphics::axis(2, at = ax2[ax2 >= 0], las = las, cex.axis=0.75)
  graphics::mtext(text = "Height", side = 2, at = mean(ax2[ax2 >= 0]), 
                  line = 3, cex = 0.95)
  
  #
  groups <- factor(levels(fac), levels = levels(fac))
  legend("topright", as.vector(groups), pch = 15, 
         col = fullColors[(1 + as.numeric(groups))], 
         bty = "n", cex = lab.cex)
  for (step in 1:nrow(mgs)) {
    for (i in 1:2) {
      x <- mgs[step, i]
      if (x < 0) {
        if(hang > 0) {
          lower.yi <- dndr_gram[step, 3] - lv_len
          label.yi <- dndr_gram[step, 3] - (1.25 * lv_len)
        } else {
          lower.yi <- min(dndr_gram[, 3]) - lv_len
          label.yi <- min(dndr_gram[, 3]) - (1.25 * lv_len)
        }
        graphics::segments(x * (-1), dndr_gram[step, 3], 
                 x * (-1), lower.yi, 
                 col = if (length(colors) ==  1) {
                   fullColors[1]
                 } else {
                   fullColors[mgs[step, 6] ]
                 }, lwd = lwd)
        graphics::text((-1) * x, label.yi,  
             #dndr_gram[step, 3] - min((1.4 * lv_len), 50), 
             hc$labels[hc$order][x * (-1)], 
             adj = c(1, 0.5), srt = 90, cex = lab.cex, 
             col = fullColors[(1 + as.numeric(fac[hc$order][x * (-1)]))])
      } else {
        x12 <- mgs[x, 3]
        y1 <- dndr_gram[step, 3]
        y2 <- dndr_gram[x, 3]
        if (length(colors) == 1) {
          colr <- fullColors[1]
        } else {
          colr <- fullColors[ mgs[step, 6] ]
        }
        graphics::segments(x12, y1, x12, y2, col = colr, lwd = lwd)
      }
    }
  }
  for (step in 1:nrow(dndr_gram)) {
    if (length(colors) == 1) {
      colr <- fullColors[1]
    } else {
      colr <- fullColors[mgs[step, 6]]
    }
    graphics::segments(dndr_gram[step, 1], dndr_gram[step, 3], 
                       dndr_gram[step, 2], dndr_gram[step, 3], 
                       col = colr, lwd = lwd)
  }
  graphics::par(mar = old.mar)
}




#' Introduction to the colorhcplot Package
#' 
#' This is a simple one-function package. 
#' Please, refer to the colorhcplot() function manual to 
#' check how the function works.
#' 
#' @details This package contains the function colorhcplot. This 
#' function generates simple colorful dendrograms and requires only 2 
#' mandatory arguments: hc and fac. The argument hc is the result of a 
#' hclust() call, while fac is a factor defining the groups. 
#' Therefore, the number of leaves of the dendrogram has to be identical to 
#' the length of fac (i.e., length(hc$labels) == length(fac) has to be TRUE). 
#' The function colorhcplot() employs a custom color palette. 
#' However, users can specify a custom list of colors.  
#' 
#' @author \packageAuthor{colorhcplot}
#' 
#' @seealso \code{\link{colorhcplot}}
#' 
#' 
#' @name colorhcplot-package
"_PACKAGE"









#' Example Gene Expression Dataset
#' 
#' This is a gene expression dataset simulating information 
#' about 499 gene probes and 13 samples, from an Affymetrix U95v2 chip. 
#' Data are made up, as well as sample labels. 
#' This dataset is adapted from the Biobase-package, version 2.32.0.
#' 
#' @format A list with 2 elements:
#' \describe{
#'   \item{exprs}{A matrix with 499 rows (genes) and 13 columns (samples)
#'     containing normalized gene expression values.}
#'   \item{fac}{A factor including the grouping for each sample.}
#' }
#' 
#' @source Data were adapted from the Biobase package version 2.32.0, and 
#' prepared by the J. Ritz Laboratory (S. Chiaretti).
#' 
#'
#' @usage data("geneData")
#'
#'
#' @examples
#' data(geneData)
#' print(geneData[[1]][1:10, 1:6])
#'
"geneData"



