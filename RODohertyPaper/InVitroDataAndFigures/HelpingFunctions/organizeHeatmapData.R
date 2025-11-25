
library(ggpubr)
library(ggplot2)
library(grid)
library(scales)
library(nleqslv)
library(tidyr)
library(exact2x2)
# takes a numeric vectors of rows of a featureCounts file, and a column integer >= 7 referring to the index of the sample of interest.
# then sums the hits to all these annotations and divides by 
# turns 0 reads into 0.5 reads
getCombinedExpression <- function(featureCountsPart, column) {
  retVal <- sum(as.numeric(featureCountsPart[,column]))/sum(as.numeric(featureCountsPart[,6]))
  if (retVal == 0) retVal = 0.5/sum(as.numeric(featureCountsPart[,6]))
  retVal
} 


# takes the sample name, the list of annotation indices (this is a list where each element corresponds to an HIV integration site, 
# and contains a vectors of the feature counts row indices corresponding to annotated genes either before or after the site, which will have their expression averaged
# Also takes the total mapped reads and GAPDH associated with the sample name.
# Returns a heatmap row (either for the pre or post HIV heatmap) which contains the sample name, then the expression of each gene (before or after HIV), and then the total mapped reads then the GAPDH
getGeneHeatmapRow <- function(sampleName, annotationIndices, featureCounts, totalMappedReads, GAPDH) {
  
  newRow <- sampleName
  
  for (i in 1:length(annotationIndices)) {
    newRow <- c(newRow, getCombinedExpression(featureCounts[annotationIndices[[i]],], which(featureCounts[1,] == sampleName)))
  }
  
  return(c(newRow, totalMappedReads, GAPDH))
}


heatmapRowOrder <- strsplit("Sample
FAST1-CRPPA
FAST8-CRPPA
FAST6-AHI1
FAST21-AHI1
FAST10-MYB
FAST13-MYB
FAST15-MYB
FAST16-MYB
FAST22-MYB
SLOW1-SOCS7
SLOW3-CENPC
SLOW5-CENPC
SLOW8-CENPC
SLOW4-C5orf24
SLOW14-OSBP2
SLOW24-GART
FAST3-HERC1
FAST17-HERC1
FAST18-HERC1
FAST4-ANKS1A
FAST5-RHOH
FAST9-SUPT3H
FAST11-TRBC2
FAST20-TRBC2
FAST23-SNRPB
SLOW6-TTLL5
SLOW7-RNF41
SLOW9-RTKN2
SLOW10-RTKN2
SLOW11-PTPRD
SLOW12-RBFOX2
SLOW15-PKN2
SLOW16-MACROD2
SLOW17-UBE2G2
SLOW18-MON2
SLOW19-LRP1B
SLOW20-BRWD1
SLOW22-NCOA7
SLOW23-CLASP2
SLOW25-ACSF3", "\n")[[1]]

heatmapGeneOrder <- strsplit("Sample  CRPPA	CRPPA	AHI1	AHI1	MYB	MYB	MYB	MYB	MYB	SOCS7	CENPC	CENPC	CENPC	C5orf24	OSBP2	GART	HERC1	HERC1	HERC1	ANKS1A	RHOH	SUPT3H	TRBC2	TRBC2	SNRPB	TTLL5	RNF41	RTKN2	RTKN2	PTPRD	RBFOX2	PKN2	MACROD2	UBE2G2	MON2	LRP1B	BRWD1	NCOA7	CLASP2	ACSF3", "[ \t]+")[[1]]
heatmapGeneOrder <- c(heatmapGeneOrder, "Mapped Transcript Levels", "GAPDH Levels")


### Averages rows of table that are given as a list of vectors of row indices that should be merged (the ... argument)
# After averaging, all but the first row in the vector is removed, and if removeCols is TRUE, the same indices of columns are removed.
averageRows <- function(inTable, removeCols = TRUE, ...) {
  toMerge <- list(...)[[1]]
  toRemove <- c()
  for (i in 1:length(toMerge)) {
    inTable[toMerge[[i]][1],] <- colMeans(inTable[toMerge[[i]],])
    toRemove <- c(toRemove, toMerge[[i]][-1])
  }
  if (removeCols) inTable <- inTable[,-toRemove]
  inTable[-toRemove,]
}



# Standardize expression by centering and scaling, we use mean and sd that counts each clone once (based on averaged heatmap), applied to the full heatmap.
# To avoid redundancy along diagonal, only show each row's gene-specific expression once (NA out the replicates). This is only done if not HIV heatmap; not relevant there.
# note we should have removed gapdh and total mapped reads columns so these are squares.
standardizeExpression <- function(regHeatmap, averagedHeatmap, notHIV = TRUE) { 
  normalizedRownames <- gsub(".*-", "",rownames(regHeatmap)) # gene name associated with 
  for (i in 1:(ncol(regHeatmap))) { 
    matchedCol <- which(colnames(averagedHeatmap) == colnames(regHeatmap)[i]) # the column of the averaged heatmap matching the non-averaged heatmap column/gene.
    regHeatmap[,i] <- (regHeatmap[,i]-mean(averagedHeatmap[,matchedCol]))/sd(averagedHeatmap[,matchedCol]) # use average/sd of averagedHeatmap version
    if (notHIV) {
      for (j in 1:nrow(regHeatmap)) {
        if (normalizedRownames[j] == colnames(regHeatmap)[i] && i != j) {
          regHeatmap[j,i] <- NA
        }
      }
    }
  }
  regHeatmap
}

library(pheatmap)
saveHeatmaps <- function(firstHeatmap, lastHeatmap, HIVHeatmap, baseFile, mainName) {
  pdf(file =paste0(baseFile, "LastExon.pdf"),   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  pheatmap(cbind(HIVHeatmap[,2], lastHeatmap),main = paste0(mainName, " Last Exons"),cluster_cols = F, scale="none", 
           cluster_rows = F,treeheight_row=0, breaks=seq(-2, 4, length.out=100))
  dev.off() 
  pdf(file =paste0(baseFile, "FirstExon.pdf"),   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  pheatmap(firstHeatmap,main = paste0(mainName, " First Exons"),cluster_cols = F, scale="none", 
           cluster_rows = F,treeheight_row=0, breaks=seq(-2, 4, length.out=100))
  dev.off() 
}
