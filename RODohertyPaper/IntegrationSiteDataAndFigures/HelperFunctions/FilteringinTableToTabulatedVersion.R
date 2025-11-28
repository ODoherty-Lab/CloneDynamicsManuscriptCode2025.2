
# Takes the 'stdTable' which is the data frame with each integration site as a row. The first five columns are, in order, the participant identifier,
# the gene symbol, the orientation ("Sense" or "Antisense"), the gene significance (i.e. a boolean value which can indicate whether it's a 
# growth-related gene, or in other cases, whether it's a highly expressed gene), and the total number of times the site was detected.
# The following columns (from column 6 to final column) represent how many times the integration site was seen at each timepoint
# So, the sum of the columns 6 to the final column should equal the value in column 5. 
# The column headers of columns 6-final column should be numeric and indicate the years on ART.


# participantIDs is a vector of strings, each being a participant ID (e.g. "CT1"). The analysis will be limited to participants in this vector.

# timeEvaluationFunc is a function which takes as an input a vector of numeric times, and returns a vector of booleans of the same length.
# It is meant to allow you to only analyze some timepoints, as retrieved from column names of columns 6-final column.
# It should TRUE at vector indices of timepoints to be included, and false for those to be excluded.
# For example, the function function(vec) {vec > 6} would only include timepoints after 6 years on ART.
# the default, defaultAllTimeEvaluation, is just function(vec) {rep(TRUE, length(vec))}. That is, it returns TRUE
# at all indices, so all timepoints are included.

# Filters the stdTable to only include specified participants and timepoints which return TRUE in the timeEvaluationFunc
# Returns output of TabulateCategories on this filtered matrix.
# This is a list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in 4th column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.

tabulateParticipant <- function(stdTable, participantIDs = NULL, timeEvaluationFunc = defaultAllTimeEvaluation) {
  
  if (length(participantIDs) != 0) {
    stdTable <- stdTable[stdTable$Participant %in% participantIDs,] # filters to only be participants.
  }
  
  timeVars <- as.numeric(colnames(stdTable)[6:ncol(stdTable)]) # gets time of sampling associated with columns 6-final column based on columns names.
  
  # a vector of columns >=6 which are associated with sampling times which "pass" the timeEvaluation functions
  relevantCols <- 5 + which(timeEvaluationFunc(timeVars))
  if (length(relevantCols) == 0) {
    print("No data within time constraints")
    return(list(NULL, NULL, NULL, NULL))
  }
  
  # addend a final column of the table to be the total number of times the clone is seen, only at timepoints passing the timeEvaluationFunc
  stdTable <- cbind(stdTable, rowSums(as.data.frame(stdTable[,relevantCols]))) 
  
  # creates output; see arguments to TabulateCategories.
  TabulateCategories(DataMatrix = stdTable, CancerRow = 4, CancerValue = TRUE, NotCancerValue = FALSE,
                     SenseRow = 3, SenseValue = "Sense", AntiValue = "Antisense", SizeRow = ncol(stdTable))
}




# takes a "DataMatrix" which is a table of clone sizes, where each row is a distinct integration site.
# specifically, it must have a column indicating whether a gene is significant (the "CancerRow") with binary values
# "CancerValue" and "NotCancerValue" indicating whether or not it's a site in a significant gene.
# must have a column indicating whether sense or antisense (the "SenseRow") with binary values "SenseValue" and "AntiValue"
# Must have a numeric column at column index "SizeRow", which contains the clone size for each integration site.
# this is only 1 column; if called by tabulateParticipant, it is the column which contains the summed clone sizes across timepoints passing the "timeEvaluationFunc".

# This list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.

TabulateCategories <- function(DataMatrix, CancerRow, CancerValue, NotCancerValue, SenseRow, SenseValue,AntiValue,SizeRow) {
  #Tabulates number of clones of each size in cancer sense. A number x at position i in ObsCPlus indicates there were x cancer sense clones of size i, treating each timepoint separately.
  CPlus <- (as.matrix(DataMatrix[DataMatrix[,SenseRow] == SenseValue & DataMatrix[,CancerRow] == CancerValue,SizeRow])) 
  CMinus <- (as.matrix(DataMatrix[DataMatrix[,SenseRow] == AntiValue & DataMatrix[,CancerRow] == CancerValue,SizeRow]))
  NPlus <- (as.matrix(DataMatrix[DataMatrix[,SenseRow] == SenseValue & DataMatrix[,CancerRow] == NotCancerValue,SizeRow]))
  NMinus <- (as.matrix(DataMatrix[DataMatrix[,SenseRow] == AntiValue & DataMatrix[,CancerRow] == NotCancerValue,SizeRow]))
  
  
  
  CPlusTab <- if(length(CPlus) == 0) NULL else tabulate(CPlus)
  CMinusTab <- if(length(CMinus) == 0) NULL else tabulate(CMinus)
  NPlusTab <- if(length(NPlus) == 0) NULL else tabulate(NPlus)
  NMinusTab <- if(length(NMinus) == 0) NULL else tabulate(NMinus)
  list(CPlusTab, CMinusTab, NPlusTab, NMinusTab)
}

# default time evaluation function which includes all timepoints.
defaultAllTimeEvaluation <- function(vec) rep(TRUE, length(vec))
