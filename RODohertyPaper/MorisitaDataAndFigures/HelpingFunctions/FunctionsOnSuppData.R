
# Takes path to TCR Simple Objects RDS file, output folder path
# Saves as Poisson-corrected version RDS in output folder with name TCRSimpleObjectsPCorrected.rds
TCRRDSToCorrectedTCRRDS <- function(pathToTCRSimpleObjectsRDS, outFolderPath, verbose = T) {
  TCRs <- readRDS(pathToTCRSimpleObjectsRDS)
  TCRsCorrected <- list()
  for (i in 1:length(TCRs)) { # shouldn't take more than minute or maybe two
    if (verbose) print(paste0("working on participant #", i))
    TCRsCorrected[[i]] <- poissonCorrectSimplifiedTCRObject(TCRs[[i]]) 
  }
  names(TCRsCorrected) <- names(TCRs)
  
  saveRDS(TCRsCorrected, paste0(outFolderPath, "/TCRSimpleObjectsPCorrected.rds"))
}

# Takes character path to supp data of proviruses (Supp Data 1), and a fold path - where to save output RDS files.
# Makes list of simple data objects, one for each participant, named by participant names. 
# Saves as ProvirusSimpleObjects.rds in given path.
# also for CT2, saves a version of simple object at same folder, "CT2WMultiples.rds", which has different row for replicate timepoints.
RDSFromSuppData1 <- function(pathToSuppData, outFolderPath, verbose = T) {
  sheets <- openxlsx::getSheetNames(pathToSuppData)
  sheets <- sheets[grep("Provirus", sheets)] # TCRs we deal with separately.
  
  participants <- unlist(strsplit(sheets, "Provirus")) # remove "Provirus" from sheet names to just get the participant IDs.
  names(participants) <- sheets
  
  provirusObjects <- list()
  
  # takes each sheet of the supp data, and turns it into a simple object in a simple object list, named by participant
  for (sheet in sheets) {
    if (verbose) print(paste0("working on sheet ", sheet))
    readSheet <- readAndCleanSuppData1(pathToSuppData, sheet = sheet)
    provirusObjects[[participants[sheet]]] <- makeSimpleObject(readSheet, colnames(readSheet))
  }
  
  CT2MultiplesSheet <- readAndCleanSuppData1(pathToSuppData, sheet = "CT2Provirus", countMults = T)
  
  saveRDS(provirusObjects, paste0(outFolderPath, "/ProvirusSimpleObjects.rds"))
  saveRDS(makeSimpleObject(CT2MultiplesSheet, colnames(CT2MultiplesSheet)), paste0(outFolderPath, "/CT2WMultiples.rds"))
}

# read in supp data 1 sheet and fix col names
# takes path to supp data 1 and sheet name, and reads it in
# making sure reading in doubles of col names didn't mess them up
# if countMults, keeps s1 and s2 separate. If not, merges them. 
readAndCleanSuppData1 <- function(pathToSuppData, sheet, countMults = FALSE) {
  readSheet <- makeNumeric(readxl::read_xlsx(pathToSuppData, sheet = sheet, skip = 1))
  colnames(readSheet) <- unlist(lapply(strsplit(colnames(readSheet), " S"), function(x){x[1]}))
  colnames(readSheet) <- round(as.numeric(colnames(readSheet)), 5) # read_xlsx apparently has bad double handling and is adding noise.
  
  if (!countMults) {
    uniqueYears <- unique(colnames(readSheet))
    toDel = c() # tracks the duplicated years to remove. I guess could also use duplicated.
    
    for (year in uniqueYears) { 
      matches <- which(colnames(readSheet) == year)
      if (length(matches) > 1) {
        toDel <- c(toDel, matches[-1])
        readSheet[,matches[1]] <- rowSums(readSheet[,matches])
      }
    }
    if (length(toDel) > 0) readSheet = readSheet[,-toDel]
  }
  
  readSheet
}

# Takes character path to supp data of TCRs (Supp Data 2), and a fold path - where to save output RDS files.
# Makes (1) list of TCR simple data objects, one for each participant, named by participant names. 
# Makes (2) list of read Excel sheets directly from supp data, for permutation test, named by participant names
# Saves (1) as TCRSimpleObjects.rds and (2) as TCRFullTables.rds in given path.
RDSFromSuppData2  <- function(pathToSuppData, outFolderPath, verbose = T) {
  
  sheets <- openxlsx::getSheetNames(pathToSuppData)
  TCRs <- list()
  TCRFullTables <- list()
  
  # takes each sheet of the supp data, and turns it into a simple object in a TCR simple object list, named by participant
  for (sheet in sheets) {
    if (verbose) print(paste0("working on sheet ", sheet))
    input = readxl::read_xlsx(pathToSuppData, sheet = sheet)
    
    # some column names were spelled out for supp data that aren't in the originally processed versions; change col names back
    colnames(input)[1] = "TCRbeta" 
    # replace -Number of Positive Reps" with "-NumPositiveReps"
    colnames(input)[grep("Number of Posi", colnames(input))] <- paste(paste0(lapply(strsplit(colnames(input)[grep("Number of Posi", colnames(input))], "-"), function(x){x[1]})), "NumPositiveReps", sep="-")
    
    TCRFullTables[[sheet]] <- input
    TCRs[[sheet]] <- SuppTCRSheetToSimpleObject(readSheet = TCRFullTables[[sheet]]) 
  } 
  
  # saves
  saveRDS(TCRFullTables, paste0(outFolderPath, "/TCRFullTables.rds"))
  saveRDS(TCRs, paste0(outFolderPath, "/TCRSimpleObjects.rds"))
  
}


# Takes character path to supp data of TCRs and name of sheet/participant to be extracted, OR an actual sheet.
# returns a simple object for the specified participant's TCRs based on the referenced TCR supp data sheet.
SuppTCRSheetToSimpleObject <- function(readSheet = NULL, pathToSuppData = NULL, sheetName = NULL) {
  
  # get the sheet
  if (is.null(readSheet)) {
    if (!is.null(pathToSuppData) && !is.null(sheetName)) readSheet <- readxl::read_xlsx(pathToSuppData, sheet = sheetName)
    else return(1)
  }
  
  readSheet2 <- readSheet[,grep("NumPositiveReps", colnames(readSheet))] # only keep sheets corresponding to number replicates positive at a timepoint
  colnames(readSheet2) <- unlist(lapply(strsplit(colnames(readSheet2), "-"), function(x){x[1]})) # renames columns to just be years.
  
  readSheetMerged <- mergeTCRRepeats(readSheet2) # combines the "B" versions of TCR sequencing by just adding together their replicates.
  makeSimpleObject(readSheetMerged, as.numeric(colnames(readSheetMerged)))
}


# takes a matrix/tibble/data frame of TCR number of positive replicates from like extracted columns of TCR supp table 4. Each column is a year name
# Adds together the columns corresponding to the same year, removes the extra columns for the repeat timepoint(s) if any present, returns data frame
mergeTCRRepeats <- function(TCRSheet) {
  
  pureYears <- unlist(strsplit(colnames(TCRSheet), "B")) # removes B if in year name
  uniqueYears <- unique(pureYears)
  toDel = c() # tracks the duplicated years to remove. I guess could also use duplicated(pureYears).
  
  TCRSheet <- makeNumeric(TCRSheet) #makes it a data frame with numeric columns. 
  
  for (year in uniqueYears) { 
    matches <- which(pureYears == year)
    if (length(matches) > 1) {
      toDel <- c(toDel, matches[-1])
      TCRSheet[,matches[1]] <- rowSums(TCRSheet[,matches])
    }
  }
  
  if (length(toDel) > 0) TCRSheet = TCRSheet[,-toDel]
  
  TCRSheet
}

# takes a data frame/matrix/tibble, makes every column numeric. Use with caution as can throw errors if not coercible.
makeNumeric <- function(mat) {
  dat <- as.data.frame(mat)
  for (i in 1:ncol(dat)) {
    dat[,i] <- as.numeric(dat[,i])
  }
  dat
}

# takes a  table where each column is a year, each row is a clone, and then a numeric vector of years of length equal to #columns of table.
# returns a simple object which is list with first element is "simpleTable" and second is "newTimes".
# simpleTable: a numeric data frame where a column is a clone (colnames are just numbered indices) and each row is a timepoint, gives clone size at time.
# newTimes: a vector of the years, corresponding to rows of simpleTable.

makeSimpleObject <- function(simpleTable, years) {
  simpleObject <- list()
  simpleObject[["simpleTable"]] <- as.data.frame(t(simpleTable))
  colnames(simpleObject$simpleTable) <- 1:ncol(simpleObject$simpleTable)
  simpleObject[["newTimes"]] <- years
  simpleObject
}
