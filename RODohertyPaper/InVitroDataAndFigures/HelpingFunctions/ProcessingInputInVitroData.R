
# takes a string which is the directory where the feature counts outputs are found.
# takes a vector of samples found within that directory whose featurecounts we want to concatenate.
# takes a two-column data frame where the first column has first part of the file names and the second column has the sample names we use in the paper.
# Takes a string extension like "exon.txt" which is concatenated in the file names. For example, "12_is_3_exon.txt" is listed as 12_is_3" in nameMatch.
# returns a single data frame with each row as an annotation, several columns of annotation metadata, and then subsequent columns being the number of hits to these annotations (via FeatureCounts) in column-indicated samples. 
# Note that the order of these output sample columns is in the same order as the samples in "sampleNames", which is important to use our code - we use it to allow our files to be a little less redundant.
concatenateProcessedFilesBySampleNames <- function(dir, sampleNames, nameMatch, extension) {
  
  fileNamesPre <- nameMatch$Sample[match(sampleNames, nameMatch$Name)] # gets vector of file names associated with sample names in same order as the sampleNames
  filesVec <- paste(fileNamesPre, extension, sep = "_") # adds on the file names like "_exon.txt" or "_gene.txt"
  filePaths <- paste(dir, filesVec, sep="/")  # paths to all featureCounts files being concatenated
  
  
  # reads an annotation file from featurecounts output (each row being annotation). 6 columns of metadata, and 1 column of read hits. Note the column names are row 1, with column names V1-V7. We make the row 1 name of column 7 as the sample name.
  readInFile <- function(filePath, sampleName) { 
    out <- read.csv(filePath, sep = "\t", header = FALSE)[-1,] # note this makes column names as first row.
    out[1,7] <- sampleName
    out
  }
  
  
  # beginning of concatenated table is the first such featureCounts file, and then we'll add columns for subsequent samples.
  outputTable <- readInFile(filePaths[1], sampleNames[1])
  if (length(filePaths) == 1) return(outputTable)
  
  # takes an existing outputTable (as described at heading of function) of concatenated featureCounts files, and adds a new read-in feature counts file (from readInFile) and adds it on. Does a quick check to make sure same annotations in first 6 columns
  concatenateNewFeatureCount <- function(existingDataFrame, newFewFeatureCount) {
    
    if(!identical(existingDataFrame[,1:6], newFewFeatureCount[,1:6])) {
      print("Not uniform input files")
      return(existingDataFrame)
    }
    
    cbind(existingDataFrame, newFewFeatureCount[,7])
  }
  
  # now concatenate on each feature count output one at a time.
  for (i in 2:length(filePaths)) {
    outputTable <- concatenateNewFeatureCount(outputTable, readInFile(filePaths[i], sampleNames[i]))
  }
  
  outputTable
}











#Takes the integration location of a provirus (HIVPosition), the orientation of a gene relative to chromosome (isPlus), and separate vectors of the gene's exon start locations and its end locations (startPositions and endPositions, which must have same length).
#Returns a list of two elements, the first is the indices of exons before HIV, and the second is indices of exons after HIV. Indices are relative to the linked startPositions and endPositions vectors.
# Only uses first exon annotated for a gene from annotations database. Note that if there's multiple annotations for a gene, the exon start and end positions will increment in one direction, and then start over when it gets to the next annotated transcript from the gene.  We use the "lastExon" numeric object below indicate the last exon associated with the first annotation.
# Considers also conventions of ordering of annotation file for opposite-orientation genes relative to chromosome.
# Leaves out an exon on each side around HIV if possible, to ignore any readthrough

getExonsBeforeAfterHIVIsoform1 <- function(startPositions, isPlus, HIVPosition, endPositions) { 
  stop <- -1 #stop is exon right before virus. Exons counted starting from 1 since R is 1-indexed language.
  lastExon <- 1
  if ((HIVPosition <= startPositions[1] && isPlus == TRUE) || (HIVPosition >= endPositions[1] && isPlus == FALSE)) { #using >= and <= because if it is equal, it is because these are RNA locations and must be splicing into it.
    stop <- 0 #provirus is before first exon
  }
  else {
    stop <- 1
    for (i in 2:length(startPositions)) {
      
      ### In the first two cases, case catches when exons jump backwards indicating end of isoform
      ### in last case, we catch if HIV has become after an exon, in which case we break.
      if (startPositions[i] < startPositions[i-1] && isPlus == TRUE) { #if gene is oriented with chromosome and numbering goes backwards
        break
      }
      if (startPositions[i] > startPositions[i-1] && isPlus == FALSE) { #if gene is oriented against chromosome and numbering jumps forwards
        break
      }
      if ((HIVPosition <= startPositions[i] && isPlus == TRUE) || (HIVPosition >= endPositions[i] && isPlus == FALSE)) {
        break #The i'th exon is after HIV.
      }
      
      stop <- i
    }
  }
  for (i in 2:length(startPositions)) {
    if (startPositions[i] < startPositions[i-1] && isPlus == TRUE) { #if gene is oriented with chromosome and numbering goes backwards
      break
    }
    if (startPositions[i] > startPositions[i-1] && isPlus == FALSE) { #if gene is oriented against chromosome and numbering jumps forwards
      break
    }
    lastExon <- i
  }
  if (lastExon == 1) {
    print("ERROR")
  }
  
  preRetValue <- c(stop, stop+1) # a length-2 vector. The first index of preRetValue will be the last index (among startPositions) of the pre-HIV exons, and its second index will be the first index (among startPositions) of the post-HIV exons.
  if (stop == 0) {
    preRetValue <- c(1,1)
    print("atStart")
    print(HIVPosition)
  }
  if(stop == length(startPositions)) {
    print("err: HIV appears to be after all exons")
    preRetValue <- c(stop,stop)
  }
  # Provide some space around HIV if possible, as there's sometimes antisense readthrough. 
  if (preRetValue[1] != 1) {
    preRetValue[1] <- preRetValue[1] - 1
  }
  if (preRetValue[2] != length(startPositions)) {
    preRetValue[2] <- preRetValue[2] + 1
  }
  list(pre = c(1:preRetValue[1]), post = c(preRetValue[2]:lastExon))
}




