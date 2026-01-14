

#Will help with text processing the year column of Bushman's files.
# YearString is value in the Bushman-provided file.
# outputs a better formatted version of it.
GetYear <- function(YearString) { 
  returnYear <- 0
  if (YearString == "CD4+ - d0") {
    returnYear <- 0
  }
  
  else {
    returnYear <- as.numeric(sub('.*y', '', YearString))
  }
  return(returnYear)
}

# Provirus is positive orientation if it's in the same orientation as the gene its in. Returns "+" or "-"
# Takes as arguments Pos, its chromosomal orientation,
# and InGeneOri, the orientation of embedding gene; if this is NA as it's intergenic, then uses nearest gene ori. 
getOri <- function(Pos, InGeneOri,NearGeneOri) { 
  PosOri <- ""
  if (grepl("-", Pos)) PosOri <- "-" # sometimes formatted weird.
  else PosOri <- "+"
  
  GeneOri <- ""
  if(is.na(InGeneOri)) GeneOri <- NearGeneOri
  else GeneOri <- InGeneOri
  
  GeneOri <- substr(GeneOri, 1, 1) # sometimes formatted weird; take first element
  
  # return + if same orientation of gene and HIV; otherwise -.
  if (PosOri == GeneOri) return("+")
  else return("-")
}



#different intergenic types depending on relative orientation/placement to nearest gene.
# takes orientation (of HIV in chromosome), 
getIntergenicType <- function(Ori, Locale, Dist) { 
  returner <- ""
  if (Locale != "intergenic") {
    returner <- Locale
  }
  else if (Ori == "+" && Dist > 0){
    returner <- "intergenic1"
  }
  else if (Ori == "+" && Dist < 0){
    returner <- "intergenic2"
  }
  else if (Ori == "-" && Dist < 0){
    returner <- "intergenic3"
  }
  else if (Ori == "-" && Dist > 0){
    returner <- "intergenic4"
  }
  returner
}


#Makes two matrices with usable data, one for each patient, synthesizing data. MatrixList[[1]] is CT1, MatrixList[[2]] is CT2
# takes as an input a specific file, BushmanCombined.csv, which contains several years of integration sites from samples labeled "CT1" and "CT2"
# note that the times sent for sequencing are intentionally scrambled, with replicate timepoints blinded, with key used at bottom to decipher.
processOriginalBushmanFile <- function(Data) {
  MatrixList <- list()
  participants <- c("CT1","CT2")
  
  #Makes two matrices with usable data, one for each patient, synthesizing data. MatrixList[[1]] is CT1, MatrixList[[2]] is CT2
  
  for (k in 1:length(participants)) {  # makes one participant's table at a time, by cycling through rows of the input Data (each being an integration site) once for each patient
    
    Patient <- participants[k] # CT1 or CT2.
    Positions <- c() # a growing vector of distinct integration sites already seen for that participant
    Years <- 0 # a growing vector of distinct timepoint labels  already seen for that participant
    
    # this is our row-growing matrix that will be MatrixList[[1]] or [[2]]. First row is the column headers.
    # there are 7 columns of metadata about the integration site (position, Bushman-annotated gene, orientation of HIV relative to gene)
    # Distance to nearest transcript start, whether Bushman labels it is cancer-related or special), 
    # and then columns for each sampling time. Under these columns we indicate the number of times a clone was observed at that timepoint
    # Each row is a clone.
    Matrix1 <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years) 
    for (i in 1:nrow(Data)) {
      if (as.character(Data[i,1]) %in% Patient) { # first column is participant
        
        # third column is a positional encoding like "chr17+4130511"
        # enter the loop below if it's an already-seen site.
        if (as.character(Data[i,3]) %in% Positions) { 
          # we add values to positions vector and rows to output matrix in sync, but output matrix has first row as header.
          # so RowIndex is the row of the table corresponding to position.
          RowIndex <- which(as.character(Data[i,3]) == Positions) + 1
          Year <- GetYear(as.character(Data[i,2])) # cleaned up timepoint value
          if (Year %in% Years) { # if we've seen this timepoint before
            Col <- which(Year == Years) + 7 # Years contains previously seen timepoints in the order of the columns, but 7 columns of metadata first
            Matrix1[RowIndex, Col] <- as.numeric(Matrix1[RowIndex, Col]) + as.numeric(Data[i,6]) # add in the number of times integration site was seen in this row of Data at this timepoint
            if (RowIndex == 1) print("d") # shouldn't happen
          }
          else { # if haven't seen timepoint before
            if (sum(Year > Years) == length(Years)) { # identical to below loop where is annotated.
              insertIndex <- length(Years)
              Matrix1 <- cbind(Matrix1, (rep(0, nrow(Matrix1))))
              Years <- c(Years, Year)
            }
            else {
              insertIndex <- max(which(Year > Years))
              Matrix1 <- cbind(Matrix1[,1:(7+insertIndex)], (rep(0, nrow(Matrix1))), Matrix1[,(insertIndex+8):ncol(Matrix1)])
              Years <- c(Years[1:insertIndex], Year, Years[(insertIndex+1):length(Years)])
              
            }
            Matrix1[1,] <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years)
            Matrix1[RowIndex, (insertIndex + 8)] <- as.numeric(Data[i,6])
            if (RowIndex == 1) print("c")
          }
          
        }
        else { # if integration site is new; hasn't been seen.
          Positions <- c(Positions, as.character(Data[i, 3])) # add integration site to the vector of positions seen.
          # generate a new row containing columns in order for the integration site (like "chr17+4130511"),
          # the Bushman-annotated gene (like "ZZEF1 *"), the relative orientation between HIV and the host gene (or nearest gene if intergenic)
          # the distance to the nearest gene transcription start (I believe this is column 8), a label of either "intronic" "intergenic" or "exonic",
          # Bushman's annotation of whether cancer-related, Bushman's annotation of whether a lymphocyte cancer-related.
          # Finally, after these 7 columns of integration site metadata, there's one column for each time of sampling
          # Times of sampling are stored in the numeric vector "Years"; they are coded.
          # We start the row by noting 0 in each year-column, and then we'll increment the value for the timepoint where we saw the sample to indicate it was seen at that time.
          NewRow <- c(as.character(Data[i, 3]), as.character(Data[i, 4]), getOri(as.character(Data[i,3]), as.character(Data[i,12]),as.character(Data[i,10])),
                      as.character(Data[i,8]), as.character(Data[i,15]), as.character(Data[i,16]), as.character(Data[i, 17]), rep(0, length(Years))) 
          NewRow[6] <- FALSE #We have our own methodology
          Matrix1 <- rbind(Matrix1, NewRow) # add on the new row into our growing matrix, since it's a distinct new position.
          RowIndex <- nrow(Matrix1) # the index of the newly added row, which is the nrow because it was just added at the end.
          Year <- GetYear(as.character(Data[i,2])) # cleans up the year labels we gave to Bushman et al. a little. Refers to the time of sampling of this newest integration site.
          
          if (Year %in% Years) { # if already saw samples from this timepoint before, and thus a column already exists for it in the table
            Col <- which(Year == Years) + 7 # The column index corresponding to the timepoint of sampling. Adding 7 as Years columns come after 7 metadata columns
            Matrix1[RowIndex, Col] <- as.numeric(Data[i,6]) # sets value to the clone size which is 6th column; remember it's newly observed site so there's no need to add the number of clones it was previously observed in
            if (RowIndex == 1) print("b") # shouldn't happen
          }
          else { # need to make new Year column
            # try to insert it at the correct ordered code position; note this isn't exactly the real order
            if (sum(Year > Years) == length(Years)) {   # if code position seems to be largest, insert new year column at end
              insertIndex <- length(Years) # index in Years directly before where we're inserting the new year column
              Matrix1 <- cbind(Matrix1, (rep(0, nrow(Matrix1)))) # add in a new Year column at end of matrix
              Years <- c(Years, Year) # add new Year to end of Years vector.
            }
            else {
              insertIndex <- max(which(Year > Years)) # index in Years directly before where we're inserting the new year column
              Matrix1 <- cbind(Matrix1[,1:(7+insertIndex)], (rep(0, nrow(Matrix1))), Matrix1[,(insertIndex+8):ncol(Matrix1)]) # insert column of 0s for new year right after insert index, adding 7 for 7 metadata columns
              Years <- c(Years[1:insertIndex], Year, Years[(insertIndex+1):length(Years)])
              
            }
            Matrix1[1,] <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years) # update first row of headers to have new Years columns in right order.
            Matrix1[RowIndex, (insertIndex + 8)] <- as.numeric(Data[i,6]) # the newly inserted year is at column index after that of insertIndex which is 7+insertIndex, so is at 8+insertIndex
            if (RowIndex == 1) print("a") # shouldn't happen
          }
        }
      }
      
    }
    
    MatrixList[[k]] <- Matrix1
  }
  
  Years <- c(2.38, 1.42, 8.19, 10.99, 10, 11.33, 12, 1.71, 10.09, 12.61, 13.74) # decode the years
  
  #Properly sets years since ART initiation, using information from Marilia.
  MatrixList[[1]][1,c(8:14)] <- Years[1:7]
  MatrixList[[2]][1,c(8:11)] <- Years[8:11]
  
  names(MatrixList) <- c("CT1", "CT2")
  MatrixList
}



GetYearNew <- function(sampleString) { # cleans up coded timepoints we sent from second Bushman file
  strsplit(sampleString, "_")[[1]][1]
}

library(rtracklayer)
# returns indices of LTR artifacts
# see email exchange with Bushman et al.; sometimes when LTR's don't quite match it gets a label of a new integration site.
# however these are generally the same as ones already seen, and chance of being a different provirus versus one already seen is very low.
# thus we flag the repeats, which are labeled with ".2" or .3, ..., at end, so we can remove them.
findArtifactLTRs <- function(posIDsWithLTRAnno) {
  LTRLabels <- unlist(lapply(strsplit(posIDsWithLTRAnno, "\\."), function(x) {x[2]})) # final character after .0, .1, or .2 or higher
  which(as.numeric(LTRLabels) > 1) # returns indieces of labels that are not .0 or .1
}


# note these lines are very similar to processOriginalBushmanFile so see this function for more detailed comments; we only comment on new things here.
# takes a very specific file.
processSecondBushmanFile <- function(newData) {
  
  newData <- newData[-findArtifactLTRs(newData$posid),] # remove the "different LTR" artifacts
  newData$posid <- unlist(lapply(strsplit(newData$posid, "\\."), function(x) {x[1]})) # remove annotation of LTR artifact type now that artifacts removed
  #newData <- newData[-which(newData$inExon | (!newData$inGene)),] #only use introns in genes. Now we do this later.
  
  
  
  participants <- c("CT1", "CT2", "CT3", "CT8", "CT9") # five participants; one matrix per participant in a list which will be returned
  MatrixListNew <- list()
  
  checkIntronic <- function(tableRow) { # takes a row of newData, and finds whether or not it's referencing an intronic gene.
    if (tableRow$inGene && !tableRow$inExon) {
      return("intronic") # for this we're just labeling if intronic or not; I know for original data we also labeled as if exonic.
    }
    return("notIntronic")
  }
  
  
  for (k in 1:length(participants)) { 
    
    Patient <- participants[k]
    Positions <- c()
    Years <- 0
    
    Matrix1 <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years)
    for (i in 1:nrow(newData)) {

      if (as.character(newData[i,3]) %in% Patient) {
        
        if (as.character(newData$posid[i]) %in% Positions) {
          RowIndex <- which(as.character(newData$posid[i]) == Positions) + 1
          Year <- GetYearNew(as.character(newData$sample[i])) # the year annotations are now found in the "sample" column.
          if (Year %in% Years) {
            Col <- which(Year == Years) + 7
            Matrix1[RowIndex, Col] <- as.numeric(Matrix1[RowIndex, Col]) + as.numeric(newData$sonicLengths[i]) # sonicLengths stores new abundance.
            if (RowIndex == 1) print("d")
          }
          else {
            if (sum(Year > Years) == length(Years)) {
              insertIndex <- length(Years)
              Matrix1 <- cbind(Matrix1, (rep(0, nrow(Matrix1))))
              Years <- c(Years, Year)
            }
            else {
              insertIndex <- max(which(Year > Years))
              Matrix1 <- cbind(Matrix1[,1:(7+insertIndex)], (rep(0, nrow(Matrix1))), Matrix1[,(insertIndex+8):ncol(Matrix1)])
              Years <- c(Years[1:insertIndex], Year, Years[(insertIndex+1):length(Years)])
              
            }
            Matrix1[1,] <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years)
            Matrix1[RowIndex, (insertIndex + 8)] <- as.numeric(newData$sonicLengths[i])
            if (RowIndex == 1) print("c")
          }
          
        }
        else {
          Positions <- c(Positions, as.character(newData$posid[i])) 
          # use gsub because sometimes Bushman gene annotations have extra symbols on them.
          # don't really care about bushman's annotations of onco and special for now. getOri functions same as before.
          # don't really care about distance to nearest gene so set it to 0.
          NewRow <- c(as.character(newData$posid[i]), gsub(",.*$", "", newData$nearestGene[i]), getOri(newData$posid[i], newData$nearestGeneStrand[i],newData$nearestGeneStrand[i]),
                      0, checkIntronic(newData[i,]), F, F, rep(0, length(Years))) 
          Matrix1 <- rbind(Matrix1, NewRow)
          RowIndex <- nrow(Matrix1)
          Year <- GetYearNew(as.character(newData$sample[i]))
          
          if (Year %in% Years) {
            Col <- which(Year == Years) + 7
            Matrix1[RowIndex, Col] <- as.numeric(newData$sonicLengths[i]) # sonicLengths stores new abundance.
            if (RowIndex == 1) print("b")
          }
          else {
            if (sum(Year > Years) == length(Years)) {
              insertIndex <- length(Years)
              Matrix1 <- cbind(Matrix1, (rep(0, nrow(Matrix1))))
              Years <- c(Years, Year)
            }
            else {
              insertIndex <- max(which(Year > Years))
              Matrix1 <- cbind(Matrix1[,1:(7+insertIndex)], (rep(0, nrow(Matrix1))), Matrix1[,(insertIndex+8):ncol(Matrix1)])
              Years <- c(Years[1:insertIndex], Year, Years[(insertIndex+1):length(Years)])
              
            }
            Matrix1[1,] <- c("Position", "Gene", "Orientation", "Distance", "Genic", "Onco", "Special", Years)
            Matrix1[RowIndex, (insertIndex + 8)] <- as.numeric(newData$sonicLengths[i])
            if (RowIndex == 1) print("a")
          }
        }
      }
      
    }
    
    MatrixListNew[[k]] <- Matrix1
  }
  

  names(MatrixListNew)[1:5] <- c("CT1", "CT2", "CT3", "CT8", "CT9")
  MatrixListNew
}




###### Some functions to help format posIDs AFTER removing the final .0 or .1

# takes character vectors of these posIDs from Bushman like "chr18-2943281"
# returns the chromosome, then a colon, and then the coordinate, which is the format used for lifting over. No orientation here
formatPosition <- function(posIDs) {
  formattedPos <- posIDs
  formattedPos[grep("-", posIDs)] <- gsub("-", ":", posIDs[grep("-", posIDs)])
  formattedPos[grep("\\+", posIDs)] <- gsub("\\+", ":", posIDs[grep("\\+", posIDs)])
  formattedPos
}

# takes character vectors of these posIDs from Bushman like "chr18-2943281"
# returns a UCSC-formatted GRanges object with first seqnames column being the chromosome and second ranges column being the coordinate, no strand indication saved here.
makeGRangeFromPosID <- function(posIDs) {
  formattedPos <- formatPosition(posIDs) # looks like chr18:2943281
  x <- strsplit(formattedPos, ":") # split it into chromosomes and coordinates
  gr <- GRanges( # makes this into a GRANGE in UCSC style, returns it.
    seqnames = unlist(lapply(x, function(x) {x[1]})),
    ranges = IRanges(as.numeric(unlist(lapply(x, function(x) {x[2]}))), width=1))
  seqlevelsStyle(gr) = "UCSC"
  gr
}


# takes character vectors of these posIDs from Bushman like "chr18-2943281"
# takes newLocs, which is a vector of new T2T coordinates. They should be of same length with matching indices.
# Replaces the coordinates in posIDs with those of T2T coordinates from newLocs, index by index, and returns this new 
# character vector of converted posIDs.
changePosidsToNewLocs <- function(posIDs, newLocs) {
  if (length(posIDs) != length(newLocs)) {
    print("ERROR: lengths unequal")
    return()
  }
  posOris <- grep("\\+", posIDs)
  for (i in posOris) {
    posIDs[i] <- paste0(strsplit(posIDs[i], "\\+")[[1]][1], "+", newLocs[i])
  }
  negOris <- grep("-", posIDs)
  for (i in negOris) {
    posIDs[i] <- paste0(strsplit(posIDs[i], "-")[[1]][1], "-", newLocs[i])
  }
  posIDs
}

###### MERGING new and old data

# takes a data matrix like the output of processSecondBushmanFile or processOriginalBushmanFile
# removes rows which aren't intronic integration sites, as indicated by 5th column being "intronic. Remember first row is headings.
keepIntronic <- function(mat) {
  toDel <- c()
  for (i in 2:nrow(mat)) {
    if (mat[i,5] != "intronic") toDel <- c(toDel, i)
  }
  if (length(toDel) > 0) mat <- mat[-toDel,]
  mat
}


# takes a data matrix like the output of processSecondBushmanFile or processOriginalBushmanFile
# removes rows which aren't intronic integration sites, as indicated by 5th column containing "notIntronic. Remember first row is headings.
keepIntronicAll <- function(mat) {
  toDel = grep("notIntronic", mat[,5])
  if (length(toDel) > 0) mat <- mat[-toDel,]
  mat
}


# takes outputs of processFirstBushmanFile and a chain to convert hg38 to T2T (from internet).
# converts output coordinates to be in T2T
# A few unconvertible intergenic sites are prepended with "hg38:" and the one unconvertible intronic site being brought to most likely T2T location but appended with "B" as in "chr7-44777049B" 
updateOldBushmanDataToT2T <- function(oldDataProcessed, conversionChain) {
  
  
  #Do one CT at a time due to anomalies. Bring CT1 and CT2 hg38 coordinates over to T2T
  x <- makeGRangeFromPosID(as.character(oldDataProcessed[[1]][-1,1]))
  x2 <- liftOver(x, conversionChain)  # brings hg38 coordinates, which is Bushman's original annotations, to the T2T coordinates of newer data coordinates
  unconvertible <- which(unlist(lapply(x2, function(x){length(x) != 1}))) # indices of oldDataProcessed[[1]][-1,1] which lack a conversion to T2T. 
  
  unconvertibleIntronic <- intersect(which(oldDataProcessed[[1]][-1, 5] == "intronic"), unconvertible) # only one of the uncovertible coordinates is intronic
  oldDataProcessed[[1]][1+unconvertibleIntronic,1] <- "chr7-44616891" # old coordinate gets brought into a repetitive region so no single conversion; this matches another clone and converts to same location correctly. I think repeat got deleted in reference update.
  #NewRow "chr7-44616883"  "OGDH *"       "-"  "10362"   "intronic"    "FALSE" "FALSE" "1"  "1"  "1"   "0"   "0"   "0"   "0" # Rest of them can just be autolabeled as t2t
  
  unconvertible <- unconvertible[-which(unconvertible == unconvertibleIntronic)] # we fixed the convertible one so remove from list of pathologies
  
  # we redo the liftover now that the OGDH position is fixed. However, we don't even try to convert the unconvertible coordinates which are at rows 1+uncovertible
  x <- makeGRangeFromPosID(as.character(oldDataProcessed[[1]][-c(1, 1+unconvertible),1])) 
  x2InsertCoordinates <- unlist(liftOver(x, conversionChain))@ranges@start # after lifting over, get a numeric vector of just the insert locations in T2T
  
  # we now fix all the coordinates in oldDataProcessed to be the T2T instead of hg38 coordinates
  #This is done in 3 steps: 
  # 1. fix the coordinates that converted correctly:
  # changePosidsToNewLocs takes the convertible (note we removed non-convertible) coordinates and replaces them with the coordinates in x2, which is a GRanges in the same order but with T2T coordinates
  oldDataProcessed[[1]][-c(1, unconvertible+1),1] <- changePosidsToNewLocs(as.character(oldDataProcessed[[1]][-c(1, unconvertible+1),1]), x2InsertCoordinates)
  
  # 2. For the OGDH clone, we're not sure it's identical site to the other OGDH clone, so append a "B" to avoid grouping them as the same clone.
  oldDataProcessed[[1]][1+unconvertibleIntronic,1] <- "chr7-44777049B"
  
  # 3. for the remaining intergenic non-convertible sites, just put the hg38 coordinate, but prepend this with "hg38:" so you can tell. We assume these have no clones in the new data.
  if (length(unconvertible) > 0) {
    oldDataProcessed[[1]][1+unconvertible,1] =  paste0("hg38:", oldDataProcessed[[1]][1+unconvertible,1])
  }
  
  
  # now do CT2, which is simpler as all its non-convertible sites are intergenic.
  y <- makeGRangeFromPosID(as.character(oldDataProcessed[[2]][-1,1]))
  y2 <- liftOver(y, conversionChain)
  unconvertible <- which(unlist(lapply(y2, function(x){length(x) != 1}))) # all intergenic
  y <- makeGRangeFromPosID(as.character(oldDataProcessed[[2]][-c(1, 1+unconvertible),1]))
  y2InsertCoordinates <- unlist(liftOver(y, conversionChain))@ranges@start
  oldDataProcessed[[2]][-c(1, 1+unconvertible),1] <- changePosidsToNewLocs(as.character(oldDataProcessed[[2]][-c(1, 1+unconvertible),1]), y2InsertCoordinates)
  oldDataProcessed[[2]][1+unconvertible,1] =  paste0("hg38:", oldDataProcessed[[2]][1+unconvertible,1])
  
  oldDataProcessed

}


# merges outputs of processFirstBushmanFile and processSecondBushmanFile respectively, after converting first output to T2T coordinates
# we treat the timepoints from the second bushman file as distinct from the first timepoint for now, and then merge them in a separate function
# removes 8th column of processSecondBushmanFile output as it's all 0s.
# for each newer integration site row from processSecondBushmanFile, tries to match its position ID
mergeTables <- function(table1, table2) {
  table2 <- table2[,-8] # 8th column of newer data is an empty timepoint; all 0 column.
  
  newTable <- cbind(table1, matrix(0, nrow(table1), ncol(table2)-7)) # adds a bunch of columns of 0s representing timepoints from newer bushman data
  for (i in 2:nrow(table2)) { # remember first row is just column titles
    rowIndex <- which(table2[i,1] == newTable[,1]) # row index the new site table2[i,1] was seen before, if it has been.
    if (length(rowIndex) == 0) { # if the new data integration site has not been seen before (converted to T2T), creates new row.
      # the new row has 7 metadata columns from table2. It then has 0s for all timepoints associated with first Bushman sequencing (=ncol(table1)-2), and then the timepoints corresponding to number of times observed in table2
      newTable <- rbind(newTable, c(table2[i, 1:7], rep(0, ncol(table1)-7), table2[i, 8:ncol(table2)])) 
    }
    else if (length(rowIndex) > 1) print(c("ERROR", i)) # shouldn't happen
    else { # if the new data integration site has been seen before (converted to T2T), combines it into old row. 
      newTable[rowIndex, (ncol(table1)+1):ncol(newTable)] <- table2[i, 8:ncol(table2)] # 1:ncol(table1 are metadata (which doesn't change) and observations in timepoints from first Bushman sampling (also unchanged). 
    }
  }
  newTable
}


# Is for a table1 like the output of processFirstBushmanFile. 
# also takes two numeric vectors cols1 and cols2 of same length.
# for each matched index of cols1 and cols2, sums the clone sizes in cols1 and cols2, and deletes the cols2.
mergeCols <- function(table1, cols, cols2) {
  if (length(cols) != length(cols2)) {
    print("ERROR")
    return()
  }
  
  for (i in 1:length(cols)) {
    table1[-1, cols[i]] <- as.numeric(table1[-1, cols[i]]) + as.numeric(table1[-1, cols2[i]])
  }
  table1 <- table1[,-cols2]
  table1
}

# merges outputs of processFirstBushmanFile and processSecondBushmanFile respectively, after converting first output to T2T coordinates
# we treat the timepoints from the second bushman file as distinct from the first timepoint for now, and then merge them in a separate function
# removes 8th column of processSecondBushmanFile output as it's all 0s.
# for each newer integration site row from processSecondBushmanFile, tries to match its position ID
# different from mergeTables because if same site is intronic in table1 and intergenic in table2, or vice versa, 
# mergeTables will label them both with the annotation of table1, whereas this will put the two sites on different rows
# with table1's clone labeled "intronic (earlier hg38-associated annotations)" and table2's clone labeled "notIntronic" (later T2T-associated annotations)"
# to help witht his, table1 "exonic"/"intergenic" get turned to "notIntronic"
mergeTablesSeparateAnnotations <- function(table1, table2) {
  table2 <- table2[,-8] # 8th column of newer data is an empty timepoint; all 0 column.
  
  newTable <- cbind(table1, matrix(0, nrow(table1), ncol(table2)-7)) # adds a bunch of columns of 0s representing timepoints from newer bushman data
  notIntronicRows <- as.numeric(which(newTable[,5] != "intronic"))[-1] # remove first row as this is just column labels
  newTable[notIntronicRows,5] <- "notIntronic"
  
  for (i in 2:nrow(table2)) { # remember first row is just column titles
    rowIndex <- which(table2[i,1] == newTable[,1]) # row index the new site table2[i,1] was seen before, if it has been.
    if (length(rowIndex) == 0) { # if the new data integration site has not been seen before (converted to T2T), creates new row.
      # the new row has 7 metadata columns from table2. It then has 0s for all timepoints associated with first Bushman sequencing (=ncol(table1)-2), and then the timepoints corresponding to number of times observed in table2
      newTable <- rbind(newTable, c(table2[i, 1:7], rep(0, ncol(table1)-7), table2[i, 8:ncol(table2)])) 
    }
    else if (length(rowIndex) > 1) print(c("ERROR", i)) # shouldn't happen
    else { # if the new data integration site has been seen before (converted to T2T), combines it into old row IF same version of intronic vs notIntronic
      if (newTable[rowIndex,5] == table2[i,5]) newTable[rowIndex, (ncol(table1)+1):ncol(newTable)] <- table2[i, 8:ncol(table2)] # 1:ncol(table1 are metadata (which doesn't change) and observations in timepoints from first Bushman sampling (also unchanged). 
      
      else {
        newTable[rowIndex, 5] <- paste0(newTable[rowIndex, 5], " (earlier hg38-associated annotations)")
        table2[i, 5] <- paste0(table2[i, 5], " (later T2T-associated annotations)")
        newTable <- rbind(newTable, c(table2[i, 1:7], rep(0, ncol(table1)-7), table2[i, 8:ncol(table2)])) 
      }
      
      
    }
  }
  newTable
}




# takes list of pMatrices like CT1CT2MergedTables or any output of Bushman processed data with the standard 7 metadata columns and then time-specific columns.
# cleans it and returns a standard matrix like in Supp Data 2 with tableToStandardParticipant (found in tableToStandardParticipant.R)

processCTMats <- function(CTMats) {
  for (i in 1:length(CTMats)) {
    CTMats[[i]][,2] <- gsub(" .*", "", CTMats[[i]][,2]) # clean gene names
    colnames(CTMats[[i]]) <- CTMats[[i]][1,] 
    CTMats[[i]] <- CTMats[[i]][-1,] # makes first row into column names
    CTMats[[i]] <- as.data.frame(CTMats[[i]]) # data frame it
    for (j in 8:ncol(CTMats[[i]])) {
      CTMats[[i]][,j] <- as.numeric(CTMats[[i]][,j]) # makes the clone size columns numeric
    }
    
    CTMats[[i]] <-  CTMats[[i]][which(CTMats[[i]][,5] == "intronic"),] #remove exon, intergenic
    CTMats[[i]][,2] <- manyAliasToSymbol(CTMats[[i]][,2]) # make the symbols standard
    
  }
  CTMats[[1]] <- cbind(rep("CT1", nrow(CTMats[[1]])), CTMats[[1]])
  CTMats[[2]] <- cbind(rep("CT2", nrow(CTMats[[2]])), CTMats[[2]]) # add participant column
  CTMats[[3]] <- cbind(rep("CT3", nrow(CTMats[[3]])), CTMats[[3]])
  CTMats[[4]] <- cbind(rep("CT8", nrow(CTMats[[4]])), CTMats[[4]])
  CTMats[[5]] <- cbind(rep("CT9", nrow(CTMats[[5]])), CTMats[[5]])
  
  ourCTsStd <- list()
  for (i in 1:5) {
    ourCTsStd[[i]] <- tableToStandardParticipant(CTMats[[i]], geneColNum = 3, oriColNum=4, sizeColNums = 9:ncol(CTMats[[i]]),
                                                 participantColNum = 1, transformAliases=FALSE)
  }
  
  mergeParticipantTables(ourCTsStd[c(1,2,3,4,5)])
}


