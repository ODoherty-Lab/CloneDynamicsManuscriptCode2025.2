
# takes a "tabulatedVersion" list. 
# This list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.


# removeUpToAndIncluding are clone sizes which are not considered large clones for calculating the clonal odds and similar values.
# for example, setting it to 1 says the size 2+ clones are used in clone odds. 2 means only size 3+ clones are considered substantial clones.
# The singleton denominator for this clone odds is always the size-1 clones regardless of this value.

# returns a table with four rows equal to the four vectors in the "tabulatedVersion" list. 
# This is like Figure 7 Table
# Also contains columns which, for each row, calculate the number of integration sites in large clones (larger than removeUpToAndIncluding)),
# The number of distinct large clones (larger than removeUpToAndIncluding), the number of clones of size 1, 
# the expansion and clone odds (defined below), and the average clone size,.
getUnformattedTable <- function(tabulatedVersion, removeUpToAndIncluding = 2) {
  
  uncollapsedVersion <- tabulatedToTable(tabulatedVersion, collapse=FALSE)
  collapsedVersion <- tabulatedToTable(tabulatedVersion, collapse=TRUE)
  
  #Make versions of functions (described below) with desired cutoff
  sumGreaterThanXClonesSpecific <- function(sizeVec) {sumGreaterThanXClones(sizeVec, removeUpToAndIncluding)}
  GetOddsOfExpansionSpecific <- function(sizeVec) {GetOddsOfExpansion(sizeVec, removeUpToAndIncluding)}
  GetOddsOfClonalitySpecific <- function(sizeVec) {GetOddsOfClonality(sizeVec, removeUpToAndIncluding)}
  sumGreaterThanXClonesUniquelySpecific <- function(sizeVec) {sumGreaterThanXClonesUniquely(sizeVec, removeUpToAndIncluding)}
  
  
  
  # Get values for table corresponding to number of distinct large clones (greater size than removeUpToAndIncluding)
  # Get values for table corresponding to total number of sampled sites in large clones (greater size than removeUpToAndIncluding)
  # Gets expansion and clone odds (these respectively are the largeUniquely and large values divided by the number of clones seen once)
  # Average clone size is total sampled integration sites divided by number of distinct clones
  # note that these are each done separately for each row, so they are vectors of four numeric values corresponding to rows 1-4.
  largeUniquely <- as.numeric(apply(uncollapsedVersion, 1, sumGreaterThanXClonesUniquelySpecific))
  large <- as.numeric(apply(uncollapsedVersion, 1, sumGreaterThanXClonesSpecific))
  expansionOdds <- as.numeric(apply(uncollapsedVersion, 1, GetOddsOfExpansionSpecific))
  cloneOdds <- as.numeric(apply(uncollapsedVersion, 1, GetOddsOfClonalitySpecific))
  averageCloneSizes <- as.numeric(apply(uncollapsedVersion, 1, GetAverageCloneSize))

  TablePre <- collapsedVersion[2:5,] # build output table with all these features. Remember 1st column of collapsed table is just the clone sizes.
  rownames(TablePre) <- c("Sense Growth", "Antisense Growth", "Sense Other", "Antisense Other")  
  outTable <- as.data.frame(cbind(TablePre, largeUniquely, large, uncollapsedVersion[,1], 
                                  expansionOdds, cloneOdds, averageCloneSizes))
  
  colnames(outTable) <- c(collapsedVersion[1,], paste(c(">", removeUpToAndIncluding, " Unique"), collapse=""),
                          paste(c("$\\geq$", removeUpToAndIncluding+1), collapse=""),
                          "<2", "Odds of Expansion", "Odds of Clone", "Average Clone Size")
  
  
  outTable
}

# takes an unformatted table output from getUnformattedTable (i.e. four columns), and rounds all values to nearest 1/100. 
roundUnformattedTable <- function(unformattedTable) {
  for (i in 1:ncol(unformattedTable)) {
    unformattedTable[,i] <- round(unformattedTable[,i], 2)
  }
  unformattedTable
}


# allows you to remove certain columns from unformmated table.
# takes an unformatted table output from getUnformattedTable (i.e. four columns; see above for specification)
# removes columns corresponding to arguments equal to FALSE. 
removeCertainUnformattedCols <- function(unformattedTable, keepUnique = TRUE, keepNonUnique = TRUE, keepCloneOdds = TRUE, keepExpOdds = TRUE,
                                         keepAveCloneSize = TRUE) {
  
  if (keepUnique == FALSE) {
    unformattedTable <- unformattedTable[,-grep("Unique", colnames(unformattedTable))]
  }
  if (keepNonUnique == FALSE) {
    potentialRemovals <- c(grep(">", colnames(unformattedTable)),grep("<", colnames(unformattedTable)),grep("geq", colnames(unformattedTable)))
    #if (length(potentialRemovals) == 2) {
    #  potentialRemovals <- potentialRemovals[2]
    #}
    unformattedTable <- unformattedTable[,-potentialRemovals]
  }
  if (keepCloneOdds == FALSE) {
    unformattedTable <- unformattedTable[,-grep("Odds of Clone", colnames(unformattedTable))]
  }
  if (keepExpOdds == FALSE) {
    unformattedTable <- unformattedTable[,-grep("Expansion", colnames(unformattedTable))]
  }
  if (keepAveCloneSize == FALSE) {
    unformattedTable <- unformattedTable[,-grep("Size", colnames(unformattedTable))]
  }
  unformattedTable
}


# Wrapper to run getUnformattedTable on the tabulatedVersion and removeUpToAndIncluding arguments; see this function's documentation.
# Then rounds its values using roundUnformattedTable and removes columns with uneeded information using removeCertainUnformattedCols
getFormattedTable <- function(tabulatedVersion, removeUpToAndIncluding = 2,
                              keepUnique = TRUE, keepNonUnique = TRUE, keepExpOdds = TRUE, keepCloneOdds = FALSE,
                              keepAveCloneSize = TRUE) {
  
  
  unformattedTable <- getUnformattedTable(tabulatedVersion, removeUpToAndIncluding=removeUpToAndIncluding)
  unformattedTable <- roundUnformattedTable(unformattedTable)
  
  
  readyToFormatTable <- removeCertainUnformattedCols(unformattedTable,keepUnique = keepUnique, keepCloneOdds = keepCloneOdds, keepNonUnique = keepNonUnique, keepExpOdds = keepExpOdds,
                                                     keepAveCloneSize = keepAveCloneSize)
  readyToFormatTable
  
}




# takes a "tabulatedVersion" list. 
# This list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.

# also takes a boolean "collapse". If TRUE, then in the output table, remove columns for which there's no clone of the given size. 
# If FALSE, have column for every clone size, even if there's no clone of that size.
# Note a further difference is the collapsed tables have an extra first row which is the clone size associated with each column.

# The output is a table where each column is a clone size
# There's four rows, each being the number of distinct clones of the row-specified type the column-associated size. 
# The first row reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second row reflects antisense clones in significant genes. The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.
# If it's a collapsed table, there's a first row before these four rows indicating the column-associated clone size; in uncollapsed tables the clone size of a column is just the column index.
tabulatedToTable <- function(tabulatedVersion, collapse) {
  
  # Makes all four lists the same length by adding 0's
  # specifically, this is saying that even if no sense integrations in significant genes of size 40 might be observed,
  # I make my vector of clone sizes stretch out to 40, noting that there's 0 sense significant clones of size 40.
  ObservedSizes <- padVecs(tabulatedVersion)
  
  
  ObsCPlus<- ObservedSizes[[1]]  # clone sizes of sense integrations into significant (e.g. growth-related) genes. The i'th element is the number of distinct size-i clones.
  ObsCMinus <- ObservedSizes[[2]] # clone sizes of sense integrations into significant (e.g. growth-related) genes. The i'th element is the number of distinct size-i clones.
  ObsNPlus <- ObservedSizes[[3]] # clone sizes of sense integrations into nonsignificant (e.g. non-growth-related) genes. The i'th element is the number of distinct size-i clones.
  ObsNMinus<- ObservedSizes[[4]] # clone sizes of sense integrations into nonsignificant (e.g. non-growth-related) genes. The i'th element is the number of distinct size-i clones.
  EmptySpots = which(ObsCPlus + ObsCMinus +  ObsNPlus + ObsNMinus == 0) # integers where there's no observed clone of this size.
  
  if (collapse && length(EmptySpots) > 0) { # if clones and we're collapsing
    #Versions with values only represented where there are actually clones
    CloneSizesCollapsed <- c(1:length(ObsCPlus))[-EmptySpots] 
    ObsCPlusCollapsed <- ObsCPlus[-EmptySpots]
    ObsCMinusCollapsed <- ObsCMinus[-EmptySpots]
    ObsNPlusCollapsed <- ObsNPlus[-EmptySpots]
    ObsNMinusCollapsed <- ObsNMinus[-EmptySpots]
    return(rbind(CloneSizesCollapsed, ObsCPlusCollapsed, ObsCMinusCollapsed, ObsNPlusCollapsed, ObsNMinusCollapsed))
  }
  else if (collapse == TRUE) {
    rbind(1:length(ObsCPlus), ObsCPlus, ObsCMinus, ObsNPlus, ObsNMinus)
  }
  else {
    return(rbind(ObsCPlus, ObsCMinus, ObsNPlus, ObsNMinus))
  }
}





# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# Returns the number of total observed integration sites of size greater than X. 
# That is, if a clone is of size i, and i>X, then it counts as i integratino sites of size greater than X.
sumGreaterThanXClones <- function(sizeVec, X = 0) {
  GreaterThanX <- 0
  if (X >= length(sizeVec)) {
    return(0)
  }
  for (i in (X+1):length(sizeVec)) {
    GreaterThanX <- GreaterThanX + i*sizeVec[i]
  }
  if (is.na(GreaterThanX)) {
    GreaterThanX <- 0
  }
  GreaterThanX
}


# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# Returns the number of distinct observed integration sites of size greater than X. 
# That is, if a clone is of size i, and i>X, then it counts as 1 integratino sites of size greater than X.
sumGreaterThanXClonesUniquely <- function(sizeVec, X = 0) {
  if (X >= length(sizeVec)) {
    return(0)
  }
  if (X != 0) {
    GreaterThanX <- sum(sizeVec[-c(1:X)])
  }
  else {
    GreaterThanX <- sum(sizeVec)
  }
  if (is.na(GreaterThanX)) {
    GreaterThanX <- 0
  }
  GreaterThanX
}

# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# the numerator is the number of total integration sites belonging to clones of size GREATER than removeUpToAndIncluding.
# the denominator is the number of size-1 clones. Returns numerator/denominator.
GetOddsOfClonality <- function(SizeVec, removeUpToAndIncluding = 2) {
  if (length(SizeVec) <= removeUpToAndIncluding) {
    return(0)
  }
  if (SizeVec[1] == 0) {
    if (sum(SizeVec[(1+removeUpToAndIncluding):length(SizeVec)]) > 0) {
      return(Inf)
    }
    else {
      return(0)
    }
  }
  sumGreaterThanXClones(SizeVec, removeUpToAndIncluding)/SizeVec[1]
}

# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# the numerator is the number of distinct clones of size GREATER than removeUpToAndIncluding.
# the denominator is the number of size-1 clones. Returns numerator/denominator.
GetOddsOfExpansion <- function(SizeVec, removeUpToAndIncluding = 2) {
  if (length(SizeVec) <= removeUpToAndIncluding) {
    return(0)
  }
  if (SizeVec[1] == 0) {
    if (sum(SizeVec[(1+removeUpToAndIncluding):length(SizeVec)]) > 0) {
      return(Inf)
    }
    else {
      return(0)
    }
  }
  singletons <- SizeVec[1]
  GreaterThanCutoff <- sum(SizeVec[(removeUpToAndIncluding+1):length(SizeVec)])
  GreaterThanCutoff/singletons
  
}



# Takes a numeric vector of tabulated clone sizes, where the i'th element is the number of distinct clones of size i.
# Returns the average clone size = total number of sampled integration sites divided by number of distinct clones.
GetAverageCloneSize <- function(SizeVec) {
  if (sum(SizeVec) == 0) {
    return(0)
  }
  sumGreaterThanXClones(SizeVec)/sum(SizeVec)
}

