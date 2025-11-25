# Produces a data frame equivalent to a sheet of supplement data 8.
# For this, uses a simple object to include basic information about factors related to comparison of Morisita between paired pairs of timepoints
# This is essentially a subset of data from makeVennSheet but in relation to matched pairs from TCR and Provirus data sets, plus comparison P values
# if TCRs (in theory object2 could be anything), object2 should be corrected tcr object and give a TCR full table.
# Labels are just for labeling features. 
getPairwiseComparisonSheet <- function(object1, object2, object2TCRMatrix = NULL, label1 = "Provirus", label2 = "TCR", resamples = 1000) {
  alignedTimepointObjects <- restrictComparisonToOverlapTable(object1, object2, object2TCRMatrix)
  object1 <- alignedTimepointObjects[[1]]
  object2 <- alignedTimepointObjects[[2]]
  object2TCRMatrix <- alignedTimepointObjects[[3]]
  
  
  times <- round(as.numeric(object1$newTimes), 2)
  
  
  outTable <- c("Timepoint (TP) 1 Years on ART", "Timepoint (TP) 2 Years on ART", "Time Elapsed", paste0(label1, " Center M-Venn"), paste0(label2, " Center M-Venn"), paste0(label1, " Center <50% Permutation Test P-Value"))
  #outTable <- c(outTable, paste0(label1, c(" Lower 2.5% When Merged and Resampled"," Upper 2.5% When Merged and Resampled")))
  #outTable <- c(outTable, paste0(label2, c(" Lower 2.5% When Merged and Resampled"," Upper 2.5% When Merged and Resampled")))
  deltaNames <- paste0(label1, c(" Center Lower 95% Delta Interval", " Center Upper 95% Delta Interval"))
  deltaNames <- c(deltaNames, paste0(label2, c(" Center Lower 95% Delta Interval", " Center Upper 95% Delta Interval")))
  outTable <- c(outTable, deltaNames, "Chi2-test Difference in Proportion n=NumClones","Z-test Difference in Proportion")
  
  
  for (i in 1:(length(times)-1)) {
    for (j in (i+1):length(times)) {
      print(j)
      newRow <- round(c(times[i], times[j], times[j] - times[i]), 1)
      newRow <- c(newRow, MorisitaDistance(object1$simpleTable[i,], object1$simpleTable[j,], TRUE, 1, FALSE, useHorn = FALSE)$b[3])
      newRow <- c(newRow, MorisitaDistance(object2$simpleTable[i,], object2$simpleTable[j,], TRUE, 1, FALSE, useHorn = FALSE)$b[3])
      #newRow <- c(newRow, getDeltaSampSize(object1, i, j), getDeltaSampSize(object2, i, j))
      newRow <- c(newRow, morisitaResamplingSignificance(object1$simpleTable[i,], object1$simpleTable[j,], simulatedPoints = resamples)[2])
      
      
      # if (is.null(object2TCRMatrix)) newRow <- c(newRow, morisitaResamplingSignificance(object2$simpleTable[i,], object2$simpleTable[j,], report95PercentileRange = TRUE))
      #  else newRow <- c(newRow, 100*TCRResamplingsTwoTimepoints(object2TCRMatrix,year1 = years[i],year2 = years[j], participantName = "X", report95PercentileRange = TRUE,makeImage = FALSE))
      
      delta1 <- getDelta95Interval(object1, i,j)
      
      delta2 <- getDelta95Interval(object2, i,j)
      
      newRow <- c(round(c(newRow, delta1, delta2), 2))
      newRow <- c(newRow, round(chiSquareTestFromObjects(object1, object2, i, j, i, j), 4))
      newRow <- c(newRow, round(getDeltaZTestFromObjects(object1, object2, i, j, i, j),4))
      outTable <- rbind(outTable, newRow)
    }
  }
  outTable
}


# takes two simple objects generally for the same person, and a full TCR table (works if this is null too though, in which case leaves it null)
# finds which of the two objects' timepoints match
# returns a list of the three arguments restricted to their matching timepoints.
restrictComparisonToOverlapTable <- function(object1, object2, object2TCRMatrix) {
  overlaps = getCommonTimeIndices(object1, object2)
  object1New <- restrictObjectToIndices(object1, overlaps[2,])
  object2New <- restrictObjectToIndices(object2, overlaps[3,])
  object2TCRMatrixNew <- restrictTCRMatrixToObject(object2TCRMatrix, object2New)
  list(object1 = object1New, object2 = object2New, object2TCRMatrix = object2TCRMatrixNew)
}


# takes two simple objects, finds where their timepoints overlap.
# returns a table with three rows. The first row are the overlapping timepoints, the second row has the index corresponding to each timepoint
# in object 1, and the third row has the index corresponding to each timepoint in object 2.
getCommonTimeIndices <- function(object1, object2) {
  times1 <- round(as.numeric(object1$newTimes), 2)
  times2 <- round(as.numeric(object2$newTimes), 2)
  overlapTimes <- times1[which(times1 %in% times2)]
  if (length(overlapTimes)==0)return(NULL)
  overlapIndices1 <- match(overlapTimes, times1)
  overlapIndices2 <- match(overlapTimes, times2)
  rbind(overlapTimes, overlapIndices1, overlapIndices2)
}

# takes a simple object and a numeric vector of timepoint indices. 
# Returns the object with the newTimes and simple Table restricted to timepoints of these indices
restrictObjectToIndices <- function(object1, indices) {
  object1$simpleTable <- object1$simpleTable[indices,]
  object1$newTimes <- object1$newTimes[indices]
  object1
}

# takes a TCR full table and simple object which may have fewer timepoints than the TCR table.
# Returns the object with the TCR full table columns restricted to timepoints present in the object.
restrictTCRMatrixToObject <- function(TCRMatrix, TCRObject) {
  if (is.null(TCRMatrix)) return(NULL)
  relevantRowTimeLabels <- TCRObject$newTimes
  TCRMatrixColLabels <- unlist(strsplit(c("Keep", unlist(lapply(strsplit(colnames(TCRMatrix[,-1]), "-"), function(x)x[1]))), "B"))
  TCRMatrix[,which(TCRMatrixColLabels %in% c("Keep", relevantRowTimeLabels))]
}


# Wrapper to run makeVennSheet on a list of simple objects, producing one data frame like in supp data 7 for each simple object
# Runs makeVennSheet on each simple object in the list of simple objects simpleObjectsArg.
# Adds output Venn sheet/supp data 7 sheet to a list of sheets originally passed to function as vennSheets; this can be empty if none yet.
# The names(vennSheets) of the added elements of the output list are equal to the concatenation of the names of the simple object with the string addendum
# So addendum should be a short (or empty) string;  these will become sheet names in excel.
# Last 3 arguments are from makeVennSheet and are passed along.
# Adds to each venn sheet an extra column which is the slope of the logit morisita line as in figure 5, with a fixed intercept.
# if doing "includeResampleP", adding a list of fullTCRTable will lead to permutations on the replicates, whereas adding none or having list index equal to 0 (e.g. with proviruses) will lead to permutation of each provirus.
makeVennSheetsList <- function(vennSheets = list(), simpleObjectsArg, addendum, includeDeltasInterval = TRUE, resamples = 10000, includeResampleP = TRUE, fullTCRTables = NULL) {
  for (i in 1:length(simpleObjectsArg)) {
    fullTCRTable <- if (is.null(fullTCRTables)) 0 else fullTCRTables[[i]]
    # wrapper for making venn sheet.
    vennSheet <- makeVennSheet(simpleObjectsArg[[i]],includeDeltasInterval = includeDeltasInterval, resamples = resamples, includeResampleP = includeResampleP, fullTCRTable = fullTCRTable)
    
    # adding slope
    distMat <- -getDistanceMatrix(simpleObjectsArg[[i]]$simpleTable, TRUE, Adjustment = 1,TRUE, useHorn=FALSE)
    times <- getLengthsAndMeans(simpleObjectsArg[[i]]$newTimes)
    slope <- lm(as.numeric(distMat) ~ 0+times$lengths)$coefficients
    vennSheet <- cbind(vennSheet, pad(c("LogitMorisitaSlope", round(as.numeric(slope), 2)), nrow(vennSheet)))
    
    # addend to passed list
    vennSheets[[paste0(names(simpleObjectsArg)[i], addendum)]] <- vennSheet
    
    print(i) # a little verbose since can be long.
  }
  vennSheets
}



# Produces a data frame equivalent to a sheet of supplement data 7.
# For this, uses a simple object to include basic information about factors going into Morisita Venn Diagram and 
# The Morisita Venn Diagram numbers themselves. Note that if there are no cross-timepoint clone pairs, we assume 2 in the Venn Diagram.
# Can also include a delta confidence interval on the Morisita center, and a permutation test p-value based on the given number of resamples.
# if doing "includeResampleP", adding a fullTCRTable will lead to permutations on the replicates, whereas having none (e.g. with proviruses) will lead to permutation of each provirus.
makeVennSheet <- function(simpleObject, includeDeltasInterval = TRUE, includeResampleP = TRUE, resamples=1000, fullTCRTable = 0) {
  cloneTable <- simpleObject$simpleTable
  simpleObject$newTimes = as.numeric(simpleObject$newTimes)
  outTable <- c("Timepoint (TP) 1 Years on ART", "Timepoint (TP) 2 Years on ART", "Time Elapsed", "TP1 Samp Size", "TP2 Samp Size", "TP1 Clone Pairs", "2X Cross-Timepoint Clone Pairs",
                "TP2 Clone Pairs", "Left Flank M-Venn", "Center M-Venn", "Right Flank M-Venn")
  if (includeResampleP) outTable <- c(outTable, "Permutation Resampled Center 5th Percentile","Permutation-test P-value")
  if (includeDeltasInterval) outTable <- c(outTable, "Center Lower 95% Delta Interval", "Center Upper 95% Delta Interval")
  for (i in 1:(nrow(cloneTable)-1)) {
    for (j in (i+1):nrow(cloneTable)) {
      print(j)
      MorisitaVenns <- MorisitaDistance(cloneTable[i,], cloneTable[j,], TRUE, 1, FALSE)
      newRow <- round(c(round(c(simpleObject$newTimes[i], simpleObject$newTimes[j], 
                                simpleObject$newTimes[j] - simpleObject$newTimes[i]), 1),sum(cloneTable[i,]), sum(cloneTable[j,]),
                        getVennPairNumbers(cloneTable[i,], cloneTable[j,]),MorisitaVenns$b[c(1,3,2)]), 2)
      if (includeResampleP) {
        if (identical(fullTCRTable,0)) newRow <- c(newRow, round(morisitaResamplingSignificance(cloneTable[i,], cloneTable[j,], simulatedPoints = resamples),3))
        else {
          TCRresamplings <- TCRResamplingsTwoTimepoints(fullTCRTable, simpleObject$newTimes[i], simpleObject$newTimes[j], maxResamples = resamples)
          newRow <- c(newRow, round(c(as.numeric(100*quantile(x=TCRresamplings$sim, probs=c(.05))), mean(TCRresamplings$sim < TCRresamplings$original)),3))
        }
      }
      if (includeDeltasInterval) { # assumes NOT horn
        deltaEstimatedSD <- sqrt(deltaMorisitaVariance(cloneTable[i,], cloneTable[j,]))
        newRow <- c(newRow, round(max(0,100*(getMorisita(cloneTable[i,], cloneTable[j,]) - 1.96*deltaEstimatedSD)),2),round(100*(getMorisita(cloneTable[i,], cloneTable[j,]) + 1.96*deltaEstimatedSD),2))
      }
      
      outTable <- rbind(outTable, newRow)
    }
  }
  outTable
}

# gives total number of sequences going into each flank of morisita venn diagram (total number of pairs of identical sequences possible in two timepoints)
# takes in two clone size vectors of equal lengths (corresponding indices), and outputs a single integer
getVennPairNumbers <- function(clones1, clones2) {
  c(sum(clones1*(clones1-1)), 2*sum(clones1*clones2), sum(clones2*(clones2-1)))
}


# helpful function that just adds "" to a vector to make it reach a certain length
pad <- function(preVec, totalLen) {
  if (length(preVec) < totalLen) {
    preVec <- c(preVec, rep("", totalLen-length(preVec)))
  }
  preVec
}
