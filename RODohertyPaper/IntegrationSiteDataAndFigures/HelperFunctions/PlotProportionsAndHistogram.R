
convertPsToStars <- function(pVec) {
  if (length(pVec) == 0) {
    return(NULL)
  }
  outVec <- rep("", length(pVec))
  outVec[pVec <= .05] <- "*"
  outVec[pVec <= .01] <- "**"
  outVec[pVec <= .001] <- "***"
  outVec
}

# takes proportions, binomials, and cutoffs.
# these are numeric vectors of the same length.
# proportions is, for each clone size cutoff, the proportion of distinct clones (generally after filtering to be in growth-related or highly-expressed genes)
# that are in the sense orientation. The i'th element considers only clones of size greater or equal to the i'th cutoff.
# binomials is a vector of p-values indicating whether the proportions of the same index are not equal to 50% by a binomial test.

# sampleName and filepath are used for labeling and naming. 
# polrOutput is the return-value of running performPolrTest on the data.
# This does a cumulative logit regression of clone size on orientation. It returns a numeric vector 3 elements, the first two being regression coefficients and the third as a p-value

# plotAntisense makes the bars represent percentAntisense instead of sense (i.e. it plots 1-proportions)
plotProportions <- function(proportions, binomials, cutoffs, sampleName = "", filepath = "", polrOutput = NULL, plotAntisense = FALSE) {
  
  if (length(binomials) != length(proportions)) {
    print("binomials and proportions must be same length")
    return(1)
  }
  
  binomialStars <- convertPsToStars(binomials) # p < .05 is 1 star, p < .01 is 2, p < .001 is 3.

  graphTitle <- paste(c(sampleName, ": Proportion sense in Effective Genes"), collapse="")
  barColor <- "#a1d999ff"
  if (plotAntisense) {
    graphTitle <- paste(c(sampleName, ": Proportion antisense in Effective Genes"), collapse="")
    proportions <- 1 - proportions
    barColor <- "#6baed6ff"
  }
  
  # make the plot
  barplot(proportions ~ cutoffs, cex.main=2, cex.lab = 1.5,cex.axis = 1.5, cex.names = 1.5, ylim=c(0,1.1), col=barColor, xlab = "Cumulative Clone Size Cutoff", ylab = "Proportion Sense", main = graphTitle)
  abline(.5,0, col="#0097a8", lty="dotted", lwd=10) #Mark proportion of 0.5.
  for (i in 1:length(cutoffs)) {
    text(x=(-.5 + i*1.2),y=proportions[i],binomialStars[i] ,adj = c(.5, .38),cex=3.5, col=ifelse(plotAntisense, "#4292c6ff", "#3a8232ff")) # add stars
  }
  if (!is.null(polrOutput)) { # adds polr coefficient (note that it's log-scale), and it's p-value, to the plot.
    text(x = 1, y = .9, labels = paste(c("polr coeff: ", format(round(polrOutput[1], 2), nsmall=2)), collapse=""))
    text(x = 1, y = .85, labels = paste(c("p = ", signif(polrOutput[3], 2)), collapse=""))
    
  }
  p <- recordPlot()
  plot.new()
  
  if (filepath == "") {
    print(p)
  }
  else {
    pdf(filepath)
    print(p)
    dev.off()
  }
  
}



# takes a "tabulatedVersion" list. 
# This list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.

# filepath and graphTitle are the save path and graph title.

# collapseTo is the largest clone size on the base of the stacked histogram; all larger clones are added to 
# this clone size. That is, if it is 5, then the stacked histogram will have columns reflecting clone sizes 1,2,3,4, and 5+.

# saves a stacked histogram of the data to the specified filepath. See description in manuscript.
makeStackedHistogram <- function(tabulatedVersion, filepath="", collapseTo = 5, graphTitle = "") {
  
  # see documentation of these functions in the "TableFormatting.R" file
  tableVersion <- tabulatedToTable(tabulatedVersion, collapse = FALSE)
  collapseToXSpecific <- function(sizeVec) {collapseToX(sizeVec, X=collapseTo)} 
  Totals <- collapseToXSpecific(colSums(tableVersion))
  Totals[Totals == 0] <- 1
  collapsedToXTable <- apply(tableVersion, 1, collapseToXSpecific)
  
  # adds larger distinct clones to the total for the max bar in the stacked histogram.
  x <-collapsedToXTable[,1]/Totals #Growth-related sense fraction
  aa <- collapseToX(tableVersion[2,], X=collapseTo)/Totals #Growth-related antisense fraction
  y <- collapseToX(tableVersion[3,], X=collapseTo)/Totals #Non-growth-related sense fraction
  z <- collapseToX(tableVersion[4,], X=collapseTo)/Totals #Non-growth-related antisense fraction
  
  
  fracsTogether <- 100*(collapsedToXTable/Totals)[,c(1,3,4,2)]
  
  # formats for number printing
  heightsCAnti <- rep(100, collapseTo)
  heightsNAnti <- rowSums(fracsTogether[,1:3])
  heightsNSense <- rowSums(fracsTogether[,1:2])
  heightsCSense <- (fracsTogether[,1])
  collapsedToXTable[collapsedToXTable == 0] <- ""
  
  
  #Plot Figure 
  a <- barplot(fracsTogether ~ c(1:collapseTo), data=fracsTogether,col=c("#E5F5E0","#A1D99B","#4292C6","#6BAED6"), border = "black", ylab = "", xlab = "", cex.names=3, cex=3, main=graphTitle)
  text(x = a, y = heightsCAnti, label = collapsedToXTable[,2], pos=1, cex = 3, col= "black")
  text(x = a, y = heightsNAnti, label = collapsedToXTable[,4], pos = 1, cex = 3, col = "black")
  text(x = a, y = heightsNSense, label = collapsedToXTable[,3], pos = 1, cex = 3, col = "black")
  text(x = a, y = heightsCSense, label = collapsedToXTable[,1], adj=c(.5,.3), cex = 3, col = "black")
  p <- recordPlot()
  plot.new()
  
  if (filepath != "") {
    pdf(file=filepath)
    print(p)
    dev.off()
  }
  
}

# uses sumGreaterThanXClonesUniquely
# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# Takes X, the desired length of the output vector
# Returns a numeric vector of length equal to X. 
# The first X-1 elements of the output are the first X-1 elements of sizeVec. 
# The X'th element of the output is the X'th element of sizeVec PLUS the sum of all remaining elements of sizeVec.
collapseToX <- function(sizeVec, X=8) {
  if (X < 1) {
    return("ERR: X must be at least 1")
  }
  if (X == 1) {
    return(sum(sizeVec))
  }
  
  outVec <- sizeVec[1:min(length(sizeVec), (X-1))]
  if (X <= length(sizeVec)) {
    outVec <- c(outVec, sumGreaterThanXClonesUniquely(sizeVec, X = X-1))
  }
  else {
    outVec <- c(outVec, rep(0, X-length(sizeVec)))
  }
  outVec
}
