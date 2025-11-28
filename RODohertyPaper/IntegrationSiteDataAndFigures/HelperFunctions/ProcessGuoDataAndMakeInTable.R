

# specifies the groupTableByClones to the Guo table; see groupTableByClones function.
HughesGroupingClonesSpecific <- function(inTable) {
  groupTableByClones(inTable, timeColumn = 9, labelColumn = 10, abundanceCol = 5)
}
# specifies the addRelativeOrientation to the Guo table; see addRelativeOrientation function.
HughesAddRelativeOrientationSpecific <- function(inTable) {
  addRelativeOrientation(inTable, 6, 4, removeUnmatchedSymbols = FALSE)
}
# does row sums of a character matrix by converting to numeric.
rowSumsCharacter <- function(characterMatrix) {
  if (ncol(characterMatrix) == 1) {
    return(as.numeric(as.character(characterMatrix[,1])))
  }
  rowSums(data.frame(apply(characterMatrix, 2, function(x) as.numeric(as.character(x)))))
}


# tableToStandardParticipant takes an input matrix with one integration site per row, which has a column for its abundance 
# at different timepoints at indices sizeColNums, gene symbols of the site in column at index geneColNum,
# participants in column at index participantColNum.
# if transformAliases is True, it will also convert all gene symbols to standard aliases, although ours already have been, and similarly can remove sites in genes without standard aliases.
# timeVars can specify the timepoints associated with sizeColNums, or can just get it from the column names if left NULL
# activeGeneColNum indicates whether gene is significant (e.g. growth-related or highly-expressed). If not indicated, will just be FALSE (can set later).
# returns table as in supp data 2 sheet 1.
tableToStandardParticipant <- function(inTable, geneColNum, oriColNum, sizeColNums, participantColNum = NULL, timeVars = NULL,
                                       activeGeneColNum = 0, transformAliases=TRUE,
                                       removeUnmatchedGenes = FALSE) {
  
  participantIDs <- (as.character(inTable[,participantColNum])) # row-aligned character vector of participant IDs.
  
  if (length(timeVars) == 0) {# if left null, sets time associated with the sizeColNums to the colnames (in numeric form)
    timeVars <- as.numeric(colnames(inTable)[sizeColNums])
    if (sum(is.na(timeVars)) > 0) {
      print("Specify time better")
    }
  }
  
  if (sum(duplicated(timeVars)) > 0 || length(timeVars) != length(sizeColNums)) {
    
    print("Process time vars better")
    return(1)
  }
  
  sizeColNums <- sizeColNums[order(timeVars)] # if times associated with columns are out of order, reorder the sizeColNums to match their proper chronology
  timeVars <- timeVars[order(timeVars)] # also order the times to be chronological. This together will enforce chronology in the output table while keeping the times matched to sizeColNums
  
  smaller <- inTable[,c(geneColNum, oriColNum, sizeColNums)] # only really need the gene column, oriColNum, and sizeColNums to get output.
  
  
  #Make sure orientations are as 1 for sense, 0 for antisense, if they're instead as "+" and "-". Remember same as TRUE/FALSE
  if ("+" %in% smaller[,2] || "-" %in% smaller[,2] ) {
    smaller[,2] <- as.numeric(smaller[,2] == "+")
  }
  smaller[,2] <- as.numeric(smaller[,2]) # make from boolean to numeric.
  
  for (i in 3:ncol(smaller)) { # for the timepoint-specific clone size columns
    smaller[,i] <- as.numeric(smaller[,i]) #make sure size columns are numeric
  }
  
  if (transformAliases == TRUE) {
    smaller[,1] <- manyAliasToSymbol(gsub(" .*", "", smaller[,1]))
  }
  
  
  
  ##significant gene labeling
  if (activeGeneColNum != 0)  smallerWithActiveGenes <- cbind(inTable[,activeGeneColNum], smaller)
  else smallerWithActiveGenes <- cbind(rep(FALSE, nrow(smaller)), smaller)
  
  # make final column order.
  smallerWithActiveGenes <- cbind(as.character(participantIDs), smallerWithActiveGenes[,2:3], smallerWithActiveGenes[,1],  
                                  smallerWithActiveGenes[,4:ncol(smallerWithActiveGenes)])
  
  colnames(smallerWithActiveGenes) <- c("Participant", "Gene", "Orientation", "sigEffect", timeVars)
  
  #if multiple times (or even if not, to be consistent), add a column combining them
  allTimes <- rowSumsCharacter(as.data.frame(smallerWithActiveGenes[,5:ncol(smallerWithActiveGenes)]))
  smallerWithActiveGenes <- as.data.frame(cbind(smallerWithActiveGenes[,1:4], allTimes, smallerWithActiveGenes[,5:ncol(smallerWithActiveGenes)]))
  for (i in 5:ncol(smallerWithActiveGenes)) {
    smallerWithActiveGenes[,i] <- as.numeric(smallerWithActiveGenes[,i])
  }
  colnames(smallerWithActiveGenes) <- c("Participant", "Gene", "Orientation", "sigEffect", "All", timeVars)
  
  smallerWithActiveGenes
  
}


# takes in a matrix where each row is an integration site, sampling time is provided in the timeColumn indexed column, 
# or else all times are assumed to be equal to timeVals. labelColumn gives the index of column with a unique clone identifier 
# abundanceCols gives index of column with the abundance sampled in the inTable row. If not given, abundances are set to abundanceVals, or 1 if not provided.
# if no label column provided, it's assumed to be last column. 

# outputs similarly structured matrix which simply groups the rows with the same integration site and sums their abundances.
groupTableByClones <- function(inTable, timeColumn = NULL, timeVals = 1, labelColumn = 0, abundanceVals = 1, abundanceCol = NULL) {
  
  if (labelColumn == 0) {
    labelColumn <- ncol(inTable)
  }
  if (length(timeColumn) == 1) {
    timeVals <- (as.numeric(inTable[,timeColumn])) # vector of times of sampling aligned with rows of table
  }
  if (length(timeVals) != nrow(inTable)) {
    timeVals <- rep(-1, nrow(inTable)) # if for whatever reason not all the timeVals are specified, just set all times to -1.
  }
  if (length(abundanceCol) == 1) {
    abundanceVals <- (as.numeric(inTable[,abundanceCol])) # vector of abundances at sampling aligned with rows of table
  }
  if (length(abundanceVals) != nrow(inTable)) {
    print("using default abundance of 1 per row")
    abundanceVals <- rep(1, nrow(inTable))
  }
  
  colsToKeep <- (1:ncol(inTable))[-c(labelColumn, timeColumn)] # indices of columns excluding the label and time columns.
  
  uniqueTimepoints <- sort(unique(timeVals)) # try to keep columns in order for timepoints. Will be new output columns
  uniqueIDs <- unique(inTable[,labelColumn]) # these will be rows of output table.
  
  
  newTable <- matrix(0, length(uniqueIDs), length(colsToKeep)+length(uniqueTimepoints)) 
  # instead of label and time columns we'll now have new column for each timepoint indicating abundance of clone at that timepoint
  
  for (i in 1:length(uniqueIDs)) { # for each uniqueID, which is a row of the output table
    relevantRows <- which(inTable[,labelColumn] == uniqueIDs[i]) # get the input table rows that reflect this clone
    newTable[i,1:length(colsToKeep)] <- as.character(inTable[relevantRows[1], colsToKeep]) # keep the metadata columns associated with this clone in the output table
    
    for (j in 1:length(relevantRows)) { # now for each "relevantRow" with the matching clone identifier, 
      abundance <- abundanceVals[relevantRows[j]] # get its abundance in the original row
      # figure out which column/timepoint is represented by the original row. This is the index of uniqueTimepoints which matches the new timepoint, plus the number of metadata columns
      colOfInterest <- which(uniqueTimepoints == timeVals[relevantRows[j]]) + length(colsToKeep) 
      newTable[i, colOfInterest] <- as.numeric(newTable[i, colOfInterest]) + as.numeric(abundance) # now increment the new table's associated timepoint-column at this clone-row by the abundance.
    }
  }
  # sums non-metadata columns for each distinct clone to get total abundance across timepoints
  totals <- rowSumsCharacter(as.data.frame(newTable[,-c(1:length(colsToKeep))])) 
  newTable <- cbind(newTable[,1:length(colsToKeep)], totals, newTable[,-c(1:length(colsToKeep))])
  
  colnames(newTable) <- c(colnames(inTable)[colsToKeep], "All", uniqueTimepoints)
  newTable
}



# takes in a matrix where each row is an integration site, and genes are present in a column and HIV strand is also in a column. Their indices are given as arguments
# 
# outputs similarly structured matrix which simply adds a final column which is the relative orientation of HIV to the host gene.
# shouldStandardizeGeneSymbols changes gene symbols in gene column to alias2Symbol standardized versions. removeUnmatchedSymbols removes rows when the gene's orientation can't be found (otherwise set to "UNK")
library(limma)
addRelativeOrientation <- function(inTable, geneCol, oriCol, shouldStandardizeGeneSymbols = TRUE, removeUnmatchedSymbols = TRUE) {
  
  txdb <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
  TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg19.knownGene
  tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")@unlistData # gene symbols/transcripts in genome
  charVers <- as.character(tx$SYMBOL) # character vector of just their symbols.
  
  
  relativeOris <- rep("", nrow(inTable)) # building new output column
  newGeneSymbols <- rep("", nrow(inTable)) # standardizing gene symbols
  
  
  for (i in 1:nrow(inTable)) { # for each row/integration site
    if (i %% 1000 == 0) print(i/nrow(inTable))
    currGene <- inTable[i,geneCol]  # get the gene name in the table
    
    if (shouldStandardizeGeneSymbols) {
      aliased <- alias2Symbol(currGene) # get its standardized version with alias2Symbol
      if (length(aliased) > 0) { # if it has a standard alias, replace it with it.
        currGene <- aliased
      }
      newGeneSymbols[i] <- currGene[1] # if multiple standard aliases use first one for updating gene symbol.
    }
    
    if (currGene[1] == "TRA") { # an unusual one which doesn't map well.
      strandOfInterest <- "+" # strandOfInterest will be the strand of the gene.
    }
    else { # if not TRA, we will find the strand of the gene based on annotations.
      hits <- c()
      for (j in 1:length(currGene)) { # go through all standard aliases from alias2Symbol
        hits <- c(hits, which(charVers == currGene[j])) # record indices of annotated genes matching this symbol.
      }
      if (length(unique(charVers[hits])) > 1) { # if there's multiple annotated symbols matching the gene symbol, print it to see how much, and use the first one
        print(unique(charVers[hits]))
        hits <- hits[1]
      }
      
      if (length(hits) == 0) { # if there's no annotated symbols matching gene symbol, set its strand to "UNK"
        print(paste(c("No matches for", currGene), collapse=""))
        strandOfInterest <- "UNK"
      }
      else {  # if there is a matching symbol (potentially after limiting to first one), get the symbol's strand annotation from UCSC genome annotations
        matchingTranscript <- tx[hits[1]]
        strandOfInterest <- as.character(matchingTranscript@strand)
      }
    }
    
    if (strandOfInterest == "UNK") { 
      relativeOris[i] <- "UNK"
    }
    else { # if know gene strand (strandOfInterest) and HIV strand (from table), set relativeOris at that index to be + if they align, and - otherwise
      relativeOris[i] <- ifelse(strandOfInterest == inTable[i, oriCol], "+", "-")
    }
    
  }
  
  
  outTable <- cbind(inTable, relativeOris) # add the new column
  if (removeUnmatchedSymbols) { # remove unknown gene symbols if desired
    toRemove <- which(relativeOris == "UNK")
    if (length(toRemove) > 0) {
      outTable <- outTable[-toRemove,]
      newGeneSymbols <- newGeneSymbols[-toRemove]
    }
  }
  
  if (shouldStandardizeGeneSymbols) { # if wanted to shouldStandardizeGeneSymbols, replace the table's gene symbol with the standardized version
    outTable[,geneCol] <- newGeneSymbols
  }
  
  outTable
}



# takes any number of tables like in supp data 2 sheet 1, and merges them by binding their rows and conglomerating timepoint columns shared between the tables
# assumes no overlapping clones between tables; different tables should be different participants for this. Otehrwise use earlier function grouping clones
mergeParticipantTables <- function(...) {
  
  tables <- list(...)[[1]]
  
  collectiveTable <- tables[[1]] # start with first table then grow it.
  
  if (length(tables) == 1) {
    print("Only 1 table input")
    
    return(collectiveTable)
  }
  
  for (i in 2:length(tables)) {
    # for each new table, figure out what time columns already exist and which new ones we'll need to makes sure we have (remember columns 6+ are timepoint specific)
    currTimeVars <- as.numeric(colnames(collectiveTable)[6:ncol(collectiveTable)]) 
    newTimeVars <- as.numeric(colnames(tables[[i]])[6:ncol(tables[[i]])])
    allTimeVars <- sort(unique(c(currTimeVars, newTimeVars))) # order them, get unique ones.
    
    priorTable <- collectiveTable # store old collectiveTable as we're going to update it.
    collectiveTable <- rbind(priorTable[,1:5], tables[[i]][,1:5]) # rbind the existing and new metadata columns of all rows for every integration site
    
    for (j in 1:length(allTimeVars)) { # now we add the time columns. Go one timepoint at a time
      # get column associated with that timepoint in the new table and existing table.
      matchingCurr <- which(currTimeVars == allTimeVars[j]) 
      matchingNew <- which(newTimeVars == allTimeVars[j])
      currPart <- rep(0, nrow(priorTable)) 
      newPart <- rep(0, nrow(tables[[i]]))
      
      if (length(matchingCurr) == 1) {# if there's a matching column in the existing data, set currPart to it. Otherwise it'll stay all 0
        currPart <- priorTable[, 5 + matchingCurr] # 5 to exclude 5 metadata columns
      }
      if (length(matchingNew) == 1) {
        newPart <- tables[[i]][, 5 + matchingNew]
      }
      
      collectiveTable <- cbind(collectiveTable, c(currPart, newPart)) # cbind on the timepoint-associated column (remember it's already time-sorted)
      # remember that the first length(currPart) rows have metaData from priorTable, while the rest have that from newer table.
    }
    
    colnames(collectiveTable) <- c("Participant", "Gene", "Orientation", "sigEffect", "All", allTimeVars) 
  }
  
  collectiveTable
}


# takes a character vector of gene symbols
# for each one, gets the alias2Symbol standard alias(es).
# if multiple hits use the first one.
# print it if it's not a pseudo numbered thing (like in some data sets) and doesn't have a hit.
# If it's "MACF1" be careful due to excel errors.
manyAliasToSymbol <- function(vec) {
  out <- rep("", length(vec))
  for (i in 1:length(vec)) {
    proposed <- alias2Symbol(vec[i])
    if (length(proposed) == 0) {
      out[i] <- vec[i]
      if(!grepl("pseudo_", vec[i])) {
        print(vec[i])
      }
    }
    else {
      out[i] <- proposed[1]
    }
    if (out[i] == "MACF1") {
      print("A mac")
      print(vec[i])
    }
  }
  out
}
