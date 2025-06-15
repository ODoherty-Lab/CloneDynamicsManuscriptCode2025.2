


##### Poisson size correction

# Poisson corrects TCR clone sizes.
# can input total number of replicates or leave it to be found (in this case leave it as NULL, in which case it defaults to 12 unless the largest clone size is <=6 in which case it's 6, or 24 if it's greater than 12).
# note maxSize (=total number replicates at a timepoint) is a vector of size equal to number of rows in the simplified clone table, i.e. number of timepoints, because each timepoint can have a different number of replicates
# returns the full TCR simple object but with the clone sizes in the simple table now adjusted
# the actual work is done by the poissonCorrectCloneSize function right below.
poissonCorrectSimplifiedTCRObject <- function(TCRSimplified, maxSize = NULL) {
  cloneTable <- TCRSimplified$simpleTable
  if (is.null(maxSize)) { # gets replicate number at each timepoint
    maxDetected <- apply(cloneTable, 1, function(x){max(as.numeric(x))})
    maxSize <- ifelse(maxDetected > 6, 12, 6)
    maxSize <- ifelse(maxDetected > 12, 24, maxSize)
    for (i in 1:length(maxDetected)) {
      if (maxDetected[i] == 3) maxSize[i] = 3
    }
  }
  
  
  for (i in 1:nrow(cloneTable)) {
    cloneTable[i,] <- poissonCorrectCloneSize(cloneTable[i,], maxSize=maxSize[i])
  }
  cloneTable <- apply(cloneTable, 2, as.numeric)
  rownames(cloneTable) <- rownames(TCRSimplified$simpleTable )
  
  TCRSimplified$simpleTable <- cloneTable
  return(TCRSimplified)
}

# takes x, a vector of TCR clone sizes that aren't adjusted, and the total replicate number maxSize
# Adjusts to method of moments poisson estimator of clone size.
# also takes maxSizeCorrection, which tells what the adjusted clone size should be for a clone seen in every replicate.
# the default values (obtained by setting to "x") are based on best-fit values assuming a power law of clone rank abundance distribution in our observed data (see Weissman et al. 2025 for more description). However findings are robust to these values.
poissonCorrectCloneSize <- function(x, maxSize = 12, maxSizeCorrection = "x") {
  if (maxSizeCorrection == "x") {
    if (maxSize == 3) maxSizeCorrection <- 6
    else if (maxSize == 6) maxSizeCorrection <- 23
    else if (maxSize == 12) maxSizeCorrection <- 54
    else if (maxSize == 24) maxSizeCorrection <- 123
  }
  x <- as.numeric(x)
  
  standardCorrections <- numeric(maxSize) 
  for (i in 1:(maxSize - 1)) {
    standardCorrections[i] <- round(-log(1-i/maxSize)*maxSize)  # Poisson method of moments estimator.
  }
  standardCorrections[maxSize] <- maxSizeCorrection
  
  # need to increment by 1 to account for size-0 clones
  standardCorrections <- c(0, standardCorrections)
  
  standardCorrections[x+1]
}


