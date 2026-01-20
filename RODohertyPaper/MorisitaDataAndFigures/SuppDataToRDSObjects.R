base = here::here()
pathToSuppData1 = paste0(base, "/MorisitaDataAndFigures/SuppData/Supplementary Data 1-ClonesTables.xlsx") # make sure copied this file into directory from paper online.
pathToSuppData2 = paste0(base, "/MorisitaDataAndFigures/SuppData/Supplementary Data 2-TCRsSequenced.xlsx") # make sure copied this file into directory from paper online.
pathToTCRRDS = paste0(base, "/MorisitaDataAndFigures/RDSObjects/TCRSimpleObjects.rds")
outFolderPath = paste0(base, "/MorisitaDataAndFigures/RDSObjects")


source(paste0(base, "/MorisitaDataAndFigures/HelpingFunctions/FunctionsOnSuppData.R")) # for supporting functions
source(paste0(base, "/MorisitaDataAndFigures/HelpingFunctions/TCRHelpFunctions.R")) # for poisson correct function

# these all should take a few min, the TCR can take few min per person because reading excel but is verbose about it, not much longer though.

# saves supp data 2 to RDS for each sheet as a list, and to list of TCR simple objects.
RDSFromSuppData2(pathToSuppData2, outFolderPath)

# saves supp data 1 to RDS for list of provirus simple objects, and for CT2 saves another RDS simple object for CT2 if timepoint repeats are counted separately
RDSFromSuppData1(pathToSuppData1, outFolderPath)


# saves conversion of RDS with list of TCR simple objects into an RDS with list of TCR simple objects where clone sizes are Poisson-corrected
TCRRDSToCorrectedTCRRDS(pathToTCRRDS, outFolderPath)
