##### Morisita Over Time Plotting Helpers

#Function that gives the pairwise distances in the order of a distance matrix... 
#kind of annoying but for a distance matrix it goes top to bottom, right to left, lower left of the diagonal.
# takes a vector of timepoints, returns a list where first element is the input vector, second element is the length of
# each pairwise interval between two times in the time vector, ordered as R would if it made a distance matrix
# and then the third element of the list is the average absolute time of each time interval with length specified in the second element
# for example, c(1,3,4) as input would be broken into ordered intervals (1,3), (1,4), and (3,4) in an R distance matrix
# so interval lengths are 2,3,1, and average times are (1+3)/2 = 2, (1+4)/2 = 2.5, and (3+4)/2 = 3.5.
getLengthsAndMeans <- function(TimeVec) {
  TimeVec = as.numeric(TimeVec)
  means <- c()
  lengths <- c()
  for (i in 1:(length(TimeVec)-1)) {
    lengths <- c(lengths, abs(TimeVec[-c(1:i)]-TimeVec[i]))
    means <- c(means, (TimeVec[-c(1:i)]+TimeVec[i])/2)
  }
  list(times=TimeVec, lengths=lengths, means=means)
}




# Graphing function which adds to a ggplot of morisita overlaps
# Takes p, a ggplot in progress.
# can add on points equal to Morisita overlap between timepoints from a simple object
# can add on trend lines of morisita overlaps regressed on time elapsed from a simple object
# to add a trendline, specify addLine = TRUE. To add points, specify addPoints = TRUE. Try to do one at a time.
# other parameters are intuitive and graphical; specify line parameters for lines and p parameters for points.
# fix0 means in line regression whether to fix intercept at morisita of 50. 
# minPoints lets it auto-recognize if you have too few points to make a good trendline and add nothing to plot.
# returns ggplot object with added feature
addToGraph <- function(p, simpleObject, distanceMat = NULL, addPoints, addLine, lineCol, pointCol = NULL, linetype,
                       adjustment = 1, minPoints = 3, psize = 4, lineSize = 4.5, pAlpha = 1, lineAlpha = 1, fix0 = TRUE) {
  

  if (length(simpleObject$simpleTable) == 0 || length(nrow(simpleObject$simpleTable)) == 0) return(p)
  if (nrow(simpleObject$simpleTable) < minPoints) return(p)
  if (is.null(distanceMat)) distanceMat <- -getDistanceMatrix(simpleObject$simpleTable, TRUE, Adjustment = adjustment,TRUE)
  
  times <- getLengthsAndMeans(simpleObject$newTimes)
  lineSlope <- lm(as.numeric(distanceMat) ~ 0+times$lengths)$coefficients
  if (!fix0) {
    lineSlope <- lm(as.numeric(distanceMat) ~ times$lengths)$coefficients[2]
  }
  
  if (addLine) {
    if (fix0) p <- p + geom_abline(intercept=0, slope = lineSlope, linetype = linetype, color=lineCol, linewidth = lineSize, alpha=lineAlpha)
    else p <- p + geom_abline(intercept=lm(as.numeric(distanceMat) ~ times$lengths)$coefficients[1], slope = lineSlope, linetype = linetype, color=lineCol, linewidth = lineSize, alpha=lineAlpha)
  }
  if (addPoints) {
    for (i in 1:length(times$lengths)) {
      p = p + geom_point(data = data.frame(cbind(times$lengths[i], distanceMat[i])),aes(x=X1,y=X2), color=pointCol, size=psize, alpha=pAlpha)
    }
  }
  p
}




###### Joint Plotting helpers

# Prints p-values from delta method and chi square method in a clean format, comparing morisita between provirus and TCR
# Specify year1 and year2 between which Morisita is being calaculated, and give two simple objects containing these timepoints. 
# figures out their year indices in each object, and calculates morisita between them in each object, and compares them.
# Prints both test p-values to console.
# Note year1b and year2b refer to timepoints in second object which can be specified different from first; if left null they're just set to year1 and year2
getMultiTestMorisitaPVals <- function(object1, object2, year1, year2, year1b = NULL, year2b = NULL) {
  if (is.null(year1b)) year1b = year1
  if (is.null(year2b)) year2b = year2
  
  i1 = which(object1$newTimes == year1)
  j1 = which(object1$newTimes == year2)
  i2 = which(object2$newTimes == year1b)
  j2 = which(object2$newTimes == year2b)
  
  if (length(i1) == 0 || length(i2) == 0 || length(j1) == 0 || length(j2) == 0) {
    print("Problem mapping years to timepoints")
    return(1)
  }
  
  zP = getDeltaZTestFromObjects(object1, object2, i1 = i1, j1 = j1, i2 = i2, j2 = j2, isObject = T)
  cP = chiSquareTestFromObjects(object1, object2, i1 = i1, j1 = j1, i2 = i2, j2 = j2, isObject = T)
  return(paste0("Chi-square: p=", signif(cP, 2), ", Delta Method: p=", signif(zP, 2)))
}



# for a proviral or TCR simple object (and full matrix if TCR), makes a row-panel with overlap venn diagram (optional), morisita venn diagram, and permutation histogram
# Input the full TCR table to make it a TCR permutation histogram. 
makeVennPanel <- function(simpleObject, year1, year2, nSims, doOverlap = T, TCRFullTable = NULL, participantName, xLims = c(0, 0.7), yLims = c(0, 4000), lwd1 = 2, lwd2 = 2, col1=col1, col2=col2, col3=col3) {
  
  row1 = which(simpleObject$newTimes == year1)
  row2 = which(simpleObject$newTimes == year2)
  
  par(mfrow = c(1,3))
  if (!doOverlap) par(mfrow = c(1,2))
  if(doOverlap) MakeVenn(simpleObject$simpleTable, row1,row2,UseNumOverlap = TRUE, col1,col2,col3, Alpha = 1) # colors defined in helper file, 6 and 7 are rows of the two 8.9 timepoints.
  MakeVenn(simpleObject$simpleTable, row1,row2,UseNumOverlap = FALSE, col1,col2,col3, Alpha = 1) # useNumOVerlap = FALSE tells it to use Morisita
  
  if (is.null(TCRFullTable)) {
    makeImageFromSharpNullProviruses(simpleObject$simpleTable[row1,],simpleObject$simpleTable[row2,], nSims = nSims, participantName, 
                                   year1, year2, xLims, yLims, lwd1, lwd2) 
  }
  
  else  makeImageFromSharpNullTCR(TCRFullTable, year1, year2, nSims, participantName, xLims, yLims, lwd1, lwd2)
  
  invisible(recordPlot())
}




###### Resampling Figures
# Takes the true Morisita between two timepoints, a vector of Morisitas from resampled or subsampled or partitioned data, participant name, the years,
#and graphical parameters.
# outputs a figure as recorded plot.
makeImageFromResampledMorisita <- function(trueMorisita, resampledMorisitas, participantName, year1, year2, xLims = c(0, 0.7), yLims = c(0, 4000), lwd1 = 2, lwd2 = 2, showAve = F) {
  hist(resampledMorisitas, main=paste0(participantName, "-", year1, year2, "Resampled"),ylim=yLims,
       xlim=xLims, breaks=seq(min(resampledMorisitas), max(resampledMorisitas), length.out=15),
       col="#16cdd0", xlab = "Permutation Morisita")
  if (!showAve) abline(v=quantile(x=resampledMorisitas, probs=c(.05)), col="black", lwd=lwd1) # black line is 5'th percentile if not showing the average
  else abline(v=mean(resampledMorisitas), col="black", lwd=lwd1) 
  abline(v=trueMorisita, col="#ff56f5ff",lwd=lwd2)
  p <- recordPlot()
  print(p)
}

# Input clone vector at two timepoints between which to get Morisita. Also get Morisita of subsamples of sizes indicates subSize1 and subSize 2 (corresponding to subsample sizes at each year)
# specify graphical parameters. Outputs a histogram figure from makeImageFromResampledMorisita showing distribution of subsamples, showing average of subsamples is near original moriista
makeImageFromSubsampledMorisitaFromData <- function(vec1, vec2, subSize1, subSize2, nSims, participantName, 
                                                    year1, year2, xLims = c(0, 0.7), yLims = c(0, 4000), lwd1 = 2, lwd2 = 2) {
  subMorisitas = getSubsampledMorisitas(vec1,
                                        vec2,
                                        subSize1, subSize2, nSims)
  trueMorisita <- getMorisita(vec1, vec2)
  makeImageFromResampledMorisita(trueMorisita, subMorisitas,  participantName, year1, year2, xLims, yLims, lwd1, lwd2, showAve = T)
}

# Input clone vector at two timepoints between which to get Morisita null distribution
# specify graphical parameters. Outputs a histogram figure from makeImageFromResampledMorisita showing distribution if you repartition between timepoints, showing 5'th percentile of null distribution in black.
makeImageFromSharpNullProviruses <- function(vec1, vec2, nSims, participantName, 
                                                    year1, year2, xLims = c(0, 0.7), yLims = c(0, 4000), lwd1 = 2, lwd2 = 2) {
  resampleMorisitas = permOneMorisitaTest(vec1, vec2, nSims)[[2]]
  trueMorisita <- getMorisita(vec1, vec2)
  makeImageFromResampledMorisita(trueMorisita, resampleMorisitas,  participantName, year1, year2, xLims, yLims, lwd1, lwd2, showAve = F)
}

# Input TCR full data from Supp Data 4, along with desired two timepoints between which to get Morisita null distribution
# specify graphical parameters. Outputs a histogram figure from makeImageFromResampledMorisita showing distribution if you repartition between timepoints, showing 5'th percentile of null distribution in black.
makeImageFromSharpNullTCR <- function(fullTCRMatrix, year1, year2, nSims, participantName, xLims = c(0, 0.7), yLims = c(0, 4000), lwd1 = 2, lwd2 = 2) {
  resampleMorisitas = TCRResamplingsTwoTimepoints(fullTCRMatrix, year1, year2, nSims)
  makeImageFromResampledMorisita(resampleMorisitas[[1]], resampleMorisitas[[2]],  participantName, year1, year2, xLims, yLims, lwd1, lwd2, showAve = F)
}



####### Venn Diagrams
col2 <- "#a0da9bb2" #right
col3 <- "#74c474b2" #center
col1 <- "#c7ebc0b2" #left
# Takes a simple table, the indices of two timepoints being compared, and the three colors of left, right, and then center flanks, and the output path (without .pdf).
# Also takes useOverlap, which tells us whether to Morisita-correct the clonal overlap (corrected when useOverlap = FALSE)
# Saves venn diagram with specified colors between specified timepoints with specified clonal overlap or morisita method to specified path.
makeAndSaveVenn <- function(simpleTable, row1, row2, useOverlap, col1, col2, col3, outPath) {
  pdf(paste0(outPath, ".pdf"))
  print(MakeVenn(simpleTable, row1, row2,useOverlap,col1, col2, col3, 1))
  dev.off()
}

library(eulerr)
library(scales)
library(grid)
library(gridBase)
# Makes Venn Digram based on a simple table (ClonesArray), the two rows/timepoints of the table being compared, whether to use NumOverlap (i.e. not Morisita method)
# also can specify color of left, right, center (in order) and a transparency factor. Can also specify whether to use Morisita-Horn instead of Morisita for Venn Diagram.
MakeVenn <- function(ClonesArray, ClonesRow1, ClonesRow2, UseNumOverlap, Col1,Col2,Col3,Alpha, useHorn = FALSE) {
  if (UseNumOverlap == TRUE) { # non-morisita method: just finds total sequences in clones that are seen at both timepoints and puts these in center
    x <- NumOverlap(ClonesArray[ClonesRow1,], ClonesArray[ClonesRow2,])
    PerA <- x$b
    PerB <- x$c
    PerBoth <- x$a
    Mid <- PerA
    Mid <- round(PerA) + round(PerBoth)
    if (Mid != 0) {
      A <- seq(1, Mid, 1)
    }
    else A <- NULL
    B <- seq((round(PerA)+1), Mid+round(PerB), 1)
  }
  else { # morisita method.
    x <- MorisitaDistance(ClonesArray[ClonesRow1,], ClonesArray[ClonesRow2,], TRUE, 1, FALSE, useHorn = useHorn)
    PerA <- x$b[1]
    PerB <- x$b[2]
    PerBoth <- x$b[3]
    
    A <- seq(1,round(PerA) + round(PerBoth), 1)
    B <- seq(round(PerA)+1, round(PerA) + round(PerBoth) + round(PerB), 1)
  }
  
  list1 <- list(NameA=A, NameB=B)
  
  names(list1) <- rownames(ClonesArray)[c(ClonesRow1,ClonesRow2)]
  
  
  plot(0, 0, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")
  
  vp <- gridBase::baseViewports()
  pushViewport(vp$figure)
  pushViewport(vp$plot)          
  on.exit(popViewport(2), add = TRUE)   
  
  
  #grid.rect(gp = gpar(fill = NA, lwd = 2))
  #plot(euler(list1), quantities=list(fontsize=18), fills = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),main = paste(c(rownames(ClonesArray)[ClonesRow1], "and",rownames(ClonesArray)[ClonesRow2]), collapse=" "))
  g <- plot(euler(list1), draw = FALSE, lwd=0, border="green",quantities=list(fontsize=18), fills = c(alpha(Col1,Alpha), alpha(Col2,Alpha), alpha(Col3,Alpha)),main = paste(c(rownames(ClonesArray)[ClonesRow1], "and",rownames(ClonesArray)[ClonesRow2]), collapse=" "))
  grid.draw(g) 
}


######## Manipulating simple objects
# takes a simple object and years to subset, returns simple object with just these years included. 
getObjectAtTimepoints <- function(simpleObject, years) {
  indices <- which(simpleObject$newTimes %in% years)
  simpleObject$simpleTable <- simpleObject$simpleTable[indices,]
  simpleObject$newTimes <- simpleObject$newTimes[indices]
  simpleObject
}

# sorts a simple table's columns by clone total clone size
sortSimpleTableCols <- function(simpleTable, cutoff = 0) {
  ordering <- order(colSums(simpleTable), decreasing = T)
  orderedTable <- simpleTable[,ordering]
  orderedTable[,colSums(orderedTable) > cutoff]
}


###### Some general R supporting functions

# takes two vectors of same length (doesn't confirm so be careful) and returns vector of this length with index-wise min between vec1 and vec2.
indexWiseMin <- function(vec1, vec2) {
  for (i in 1:length(vec1)) {
    vec1[i] <- min(vec1[i], vec2[i])
  }
  vec1
}


# runs Barnard's barnad.test function and returns the p value but doesn't print anything to console verbosely.
library(Barnard)
quiet_barnard <- function(...) {
  invisible(
    capture.output(
      result <- Barnard::barnard.test(...),           
      file = tempfile()                              
    )
  )
  result$p.value[2]                                           
}
