colorhcplot <-
function (hc, 
                         fac, 
                         hang = 0.1, 
                         main = "Cluster Dendrogram", 
                         colors = NULL, 
                         lab.cex = 1, 
                         ylim = NULL, 
                         lwd = 3,
                         las = 1,
                         lab.mar = 0.55)
{
  
  if (class(hc) != "hclust")
    stop("Please, provide an 'hclust'-class object")
  
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
  old.mar <- par()$mar
  par(mar = c(2.1, 4.1, 4.1, 2.1))
  
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
  
  ax2 <- axis(2, pos = -10000)
  axis(2, at = ax2[ax2 >= 0], las = las, cex.axis=0.75)
  mtext(text = "Height", side = 2, at = mean(ax2[ax2 >= 0]), line = 3, cex = 0.95)

  #
  groups <- factor(levels(fac), levels = levels(fac))
  legend("topright", as.vector(groups), pch = 15, col = fullColors[(1 + as.numeric(groups))], 
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
        segments(x * (-1), dndr_gram[step, 3], 
                 x * (-1), lower.yi, 
                 col = if (length(colors) ==  1) {
                   fullColors[1]
                 } else {
                   fullColors[mgs[step, 6] ]
                 }, lwd = lwd)
        text((-1) * x, label.yi,  
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
        segments(x12, y1, x12, y2, col = colr, lwd = lwd)
      }
    }
  }
  for (step in 1:nrow(dndr_gram)) {
    if (length(colors) == 1) {
      colr <- fullColors[1]
    } else {
      colr <- fullColors[mgs[step, 6]]
    }
    segments(dndr_gram[step, 1], dndr_gram[step, 3], dndr_gram[step, 
                                                               2], dndr_gram[step, 3], col = colr, lwd = lwd)
  }
  par(mar = old.mar)
}
