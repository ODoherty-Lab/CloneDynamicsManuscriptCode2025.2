# Please only load this script along with the other helper functions (which is done in the "ProduceIntegrationSiteOutputs.Rmd"
# as it references other R functions which should be loaded from there; hard to coordinate filepaths in each script).

# Takes the 'inTable' which is the data frame with each integration site as a row. The first five columns are, in order, the participant identifier,
# the gene symbol, the orientation ("Sense" or "Antisense"), the gene significance (i.e. a boolean value which can indicate whether it's a 
# growth-related gene, or in other cases, whether it's a highly expressed gene), and the total number of times the site was detected.
# The following columns (from column 6 to final column) represent how many times the integration site was seen at each timepoint
# So, the sum of the columns 6 to the final column should equal the value in column 5. 
# The column headers of columns 6-final column should be numeric and indicate the years on ART.

# The sample name is used to name the file outputs. These file outputs can also be manually provided. (the "FileName" arguments)
# autoNameWithDirectory must be given to auto-name the files as it specifies the route folder at which things will be saved.

# participantIDs is a vector of strings, each being a participant ID (e.g. "CT1"). The analysis will be limited to participants in this vector.

# timeEvaluationFunc is a function which takes as an input a vector of numeric times, and returns a vector of booleans of the same length.
# It is meant to allow you to only analyze some timepoints, as retrieved from column names of columns 6-final column.
# It should TRUE at vector indices of timepoints to be included, and false for those to be excluded.
# For example, the function function(vec) {vec > 6} would only include timepoints after 6 years on ART.
# the default, defaultAllTimeEvaluation, is just function(vec) {rep(TRUE, length(vec))}. That is, it returns TRUE
# at all indices, so all timepoints are included.

# removeUpToAndIncluding are clone sizes which are not considered large clones for calculating the clonal odds and similar values.
# for example, setting it to 1 says the size 2+ clones are used in clone odds. 2 means only size 3+ clones are considered substantial clones.
# The singleton denominator for this clone odds is always the size-1 clones regardless of this value.

# stackedHistogramCollapseTo is the largest clone size on the base of the stacked histogram; all larger clones are added to 
# this clone size. That is, if it is 5, then the stacked histogram will have columns reflecting clone sizes 1,2,3,4, and 5+.

# cutoffs gives the clone size cutoffs on the columns of the proportion sense plots.
# cutoffs is a numeric vector of any length > 1. Each value is the inclusive cutoff of one of the bars on the proportion plot.
# e.g. if cutoffs are c(1,2,4,7), you'll get bars for clone size >=1, size >=2, >=4, and >= 7.
# oneSided can be 0, 1, or 2, reflecting 2-sided, greater alternative, and lesser alternative respectively, for performing proportion plot binomial tests.

# plotAntisense plots percent antisense instead of percent sense; this isn't used in paper.
# the "keep" arguments are about columns used in the outputted clone size table; set to TRUE to include these values.
# specifically, "Unique" is the number of sites of size 1, "NonUnique" is number larger clones, ExpOdds is the 
# the number of clones divided by number of singles, counting each clone once. Clone odds is the same except it counts multiple times.
# Average clone size is the total number of sequenced sites divided by number unique sites.

# colsToSplit allows wrapping of the clone size table onto multiple rows. See documentation of outputNiceTable function but recommended to set to "auto".

# produces figures at the filepaths as a bar graph of proportions sense (or antisense if specified) with polr output
# and binmomial significance shown, a clear PDF'd Latex clone size table and also saves this as a CSV, and a "stacked histogram" plot
# Also returns a list with the 4 elements
# the first is the "tabulatedVersion" of the data used after filtering timepoints and participants by the arguments.
# this is itself a list of 4 elements described in the "tabulateParticipant" documentation
# the second is the polr output, which gives the cumulative logit regression coefficient and p-value.
# The third is the formatted table, which is data frame form of the clone size table separately saved as a CSV and PDF.
# the final is the proportionSenseOutput
# this has the proportions sense at each cutoff and its binomial significance as reflected in the proportions plot.

doEverything <- function(inTable,  sampleName, proportionFileName = "", 
                         participantIDs = NULL, timeEvaluationFunc = defaultAllTimeEvaluation,
                         cutoffs = NULL,
                         removeUpToAndIncluding = 1, 
                         stackedHistogramFileName = "",stackedHistogramCollapseTo=5,
                         formattedTableCaption = "", formattedTableFile = NULL, colsToSplit = NULL,
                         keepUnique = TRUE, keepNonUnique = TRUE, keepExpOdds = FALSE,
                         keepAveCloneSize = TRUE, keepCloneOdds = TRUE,
                         autoNameWithDirectory = NULL, plotAntisense = FALSE, oneSided = 0) {
  
  tabulatedVersion <- tabulateParticipant(inTable, participantIDs = participantIDs, timeEvaluationFunc = timeEvaluationFunc)
  
  
  if (length(autoNameWithDirectory) == 1) { # set file names/paths if autoNameWithDirectory is given as a folderpath.
    setwd(autoNameWithDirectory)
    proportionFileName <- paste(c(sampleName, "Proportions.pdf"), collapse="")
    stackedHistogramFileName <- paste(c(sampleName, "StackedHistogram.pdf"), collapse="")
    formattedTableFile <- paste(c(sampleName, "Table.pdf"), collapse="")
    formattedTableCaption <- paste(c(sampleName, "Clone Sizes"), collapse=" ")
  }
  
  ### Get proportions sense and binomial p-values
  proportionSenseOutput <- getProportionSenseAtCutoffsWithBinomials(tabulatedVersion, cutoffs = cutoffs, oneSided = oneSided)
  
  ### Perform polr test
  polrOutput <- performPolrTest(tabulatedVersion)
  
  ### MakeOutputTable
  unformattedTable <- getUnformattedTable(tabulatedVersion, removeUpToAndIncluding)
  formattedTable <- getFormattedTable(tabulatedVersion, removeUpToAndIncluding = removeUpToAndIncluding,
                                      keepUnique = keepUnique, keepCloneOdds = keepCloneOdds, keepNonUnique = keepNonUnique, keepExpOdds = keepExpOdds,
                                      keepAveCloneSize = keepAveCloneSize) 
  
  ### If saving a formatted table is desired
  if (length(formattedTableFile) == 1) {
    
    # outputNiceTable not currently implemented as I need to add comments to these tedious functions; just save the CSV version for now.
    # ordinarily formats a Latex version with nice spacing.
    #outputNiceTable(formattedTable, filepathTex = gsub(".pdf$", ".tex", formattedTableFile), caption = formattedTableCaption, colsToSplit = NULL, outputPath = autoNameWithDirectory) 
    write.csv(formattedTable, gsub(".pdf$", ".csv", formattedTableFile))
    
  }
  
  ### Plot the proportions
  plotProportions(proportionSenseOutput$proportions, proportionSenseOutput$binomials, 
                  proportionSenseOutput$cutoffsUsed, sampleName = sampleName, filepath = proportionFileName, polrOutput = polrOutput, plotAntisense = plotAntisense) 
  
  ### Make the stacked histogram
  
  makeStackedHistogram(tabulatedVersion, filepath=stackedHistogramFileName, collapseTo = stackedHistogramCollapseTo, graphTitle = sampleName) 
  
  list(tabulatedVersion, polrOutput, formattedTable, proportionSenseOutput)
  
}





# this is a wrapper that runs doEverything for different vectors of participants and autocreates output folders; refer to documentation above
# the "..." represents a list of vectors of participant IDs. 
# doEverything will be run once (using inTable = bigTable as described above) for each vector of participantIDs.
# allows specification of the folder to save everything in (directoryToUse) which will be set as autoNameWithDirectory in doEverything.
# allows specification of a limited number of the doEverything arguments (plotAntisense and oneSided); ideally the ... would allow arbitrary arguments specification but not presently.
# cutoffs are fixed to c(1,2,4,6,8), with other arguments about what to keep as specified in the doEverything call below.

# output path folder is created if it doesn't exist. 
# finally, in addition to creating figures, it returns a list with one element for each participantIDs vector.
# Each of these lists' elements represents the return value of getProportionSenseAtCutoffsWithBinomials;
# this has the proportions sense at each cutoff and its binomial significance as reflected in the proportions plot.

doEverythingForIDSet <- function(bigTable, directoryToUse, plotAntisense = FALSE,oneSided = 1, ...) { 
  participantSets <- list(...)[[1]]
  outputRegressions <- list()
  outputCIs <- list()
  dir.create(directoryToUse, showWarnings = FALSE)
  for (i in 1:length(participantSets)) {
    
    everythingOut <- doEverything(bigTable,  participantSets[[i]], autoNameWithDirectory = directoryToUse, participantIDs = participantSets[[i]], plotAntisense = plotAntisense,
                                  cutoffs=c(1,2,4,6,8), colsToSplit = "auto", keepUnique = FALSE, oneSided = oneSided, keepCloneOdds = TRUE, keepNonUnique = FALSE)
    outputCIs[[i]] <- everythingOut[[4]]
  }
  outputCIs
}
