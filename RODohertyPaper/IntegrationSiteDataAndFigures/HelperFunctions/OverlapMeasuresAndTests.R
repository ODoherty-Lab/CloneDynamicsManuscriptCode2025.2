# These functions calculate the Simpson and Morisita 

#In general, the Morisita numerator is the probability of seeing identical clones if a provirus/clone is selected from timepoint 1 and another is selected from timepoint 2.
#The denominator is the sum of the numerator with the simpson indices of timepoints 1 and 2 (where the simpson is the probability of seeing identical clones if two proviruses/TCRs are drawn from a single timpeoint).

#We also have a function here which does fast subsampling without replacement from a vector of clone sizes.

#Finally, we include an older version of NumOverlap and Morisita Distance functions that gives each part of the venn diagram pre and post Morisita which is used in the plotting functions


############# MORISITA HELPERS

### clones1 is a numeric vector of clone sizes at timepoint 1, can contain 0s.
### clones2 is a numeric vector of clone sizes at timepoint 2, can contain 0s.
### Equal indices in clones1 and clones2 should correspond to the same clone, so vectors must be equal length

getMorisitaNumerator <- function(clones1, clones2) { 
  2*sum(clones1*clones2)/(as.numeric(sum(clones1))*as.numeric(sum(clones2)))
}

Simpson2 <- function(clones1) {
  sum(clones1*(clones1-1)/(sum(clones1)^2-sum(clones1)))
}

getMorisitaDenominator <- function(clones1, clones2) { 
  simp1 <- Simpson2(clones1)
  simp2 <- Simpson2(clones2)
  simpCenter <- getMorisitaNumerator(clones1, clones2)
  simp1+simp2+simpCenter
}

# if avoid 0, assumes the morisita numerator gets no lower than a single overlapping clone.
getMorisita <- function(clones1, clones2, avoid0 = F) {
  simp1 <- Simpson2(clones1)
  simp2 <- Simpson2(clones2)
  simpCenter <- 2*sum(clones1*clones2)/(as.numeric(sum(clones1))*as.numeric(sum(clones2)))
  if (avoid0 && simpCenter == 0) simpCenter = 2/(sum(clones1)*sum(clones2))
  simpCenter/getMorisitaDenominator(clones1, clones2)
}



# takes a numeric vector of clones sizes, and quickly expands it into a large vector where each index is a TCR/provirus with a number corresponding to its clone identifier. 
# In the output, indices with identical numbers correspond to clones
# This is important for reasonable-speed subsampling without replacement. 
elongate <- function(clones) {
  clones = as.numeric(clones)
  elongated <- rep(0, sum(clones))
  max <- length(clones)
  currIndex <- 1
  for (i in 1:max) {
    if (clones[i] == 0) next
    elongated[currIndex:(currIndex+clones[i]-1)] <- i
    currIndex <- currIndex + clones[i]
  }
  elongated
}



# takes clones1 and clones2 as described earlier, 
# takes merged which should be clones1+clones2, and represents the "merged" timepoint
# takes elongated, which is just elongate(merged), in case running this many times and don't want to wait for elongate to run a bunch
# if merged or elongated aren't provided, they're calculated.
# Then, function will re-partition the merged timepoint into a new clones1 and clones2 (resampled1 and resampled2) of the same lengths as clones2 but random which clones are there
# returns the morisita between resampled1 and resampled2
# IMPORTANT: Equivalent to randomly permuting timepoint labels between proviruses/TCRs in clones1 and clones2, and then calculating morisita between the new "timepoints."
morisitaResampleWithElongated <- function(clones1, clones2, merged=NULL, elongated = NULL) {
  if (is.null(merged) || is.null(elongated)) {
    merged = clones1+clones2
    elongated <- elongate(merged)
  }
  # THIS IS CRUCIAL STEP TO QUICKLY SUBSAMPLE WITHOUT REPLACEMENT FROM THE ELONGATED VECTOR:
  resample1 <- tabulate(sample(elongated, sum(clones1), FALSE)) 
  if (length(merged) > length(resample1)) resample1 <- c(resample1, rep(0, length(merged)-length(resample1))) # pads clone size vector.
  resample2 <- merged-resample1 # the remaining clones.
  getMorisita(resample1, resample2)
}

# takes a numeric vector of clone sizes clones1, and "elongated" version of it for fast subsampling or can calculate this.
#Returns a numeric vector of clone sizes in subsample, with indices corresponding to the same clones as in clones1.
resampleFromElongated <- function(clones1, sampSize, elongated = NULL) {
  if (is.null(elongated)) elongated <- elongate(clones1)
  
  resample1 <- tabulate(sample(elongated, sampSize, FALSE)) 
  if (length(clones1) > length(resample1)) resample1 <- c(resample1, rep(0, length(clones1)-length(resample1))) # pad resample to be same length as clones1
  
  return(resample1)
}

#### For getting Morisita directly from matrix of replicate presence-absence data.
# Takes a full matrix of TCR presence-absence data 
# Here each column is a replicate, and each cell is whether or not the row-indicated TCR is present (1) or absent (0)
# Takes two numeric vectors of columns representing columns of replicates of timepoints we take the Morisita between.
# Poisson-corrects and returns Morisita between the two sets of timepoints.
morisitaFromReadCols <- function(mat, colSet1, colSet2) {
  clones1 <- poissonCorrectCloneSize(rowSums(mat[,colSet1]), maxSize = length(colSet1), maxSizeCorrection = "x")
  clones2 <- poissonCorrectCloneSize(rowSums(mat[,colSet2]), maxSize = length(colSet2), maxSizeCorrection = "x")
  getMorisita(clones1, clones2)
}

# takes paired numeric vectors of clone sizes clones1 and clones2
# takes sizes of subsamples of clones1 and clones2 desired
# subsamples of these sizes are drawn without repalcement from clones1 and clones2 a number of times equal to nSims, and Morisitas are calculated between the subsamples.
# returns vector of length nSims containing moristas between subsamples.
getSubsampledMorisitas <- function(clones1, clones2, sampSize1, sampSize2, nSims) {
  if (sum(clones1) < sampSize1 | sum(clones2) < sampSize2) {
    print("samp sizes more than total clones")
    return()
  }
  
  elongated1 = elongate(clones1)
  elongated2 = elongate(clones2)
  
  morisitas <- numeric(nSims)
  for (i in 1:nSims) {
    morisitas[i] = getMorisita(resampleFromElongated(clones1, sampSize1, elongated = elongated1), 
                               resampleFromElongated(clones2, sampSize2, elongated = elongated2))
  }
  morisitas
}


#Older function to calculate Morisita Distance. AsProportion uses the way we do it, if FALSE it's the classic Morisita dist (=1-Morisita index)
#Vec1 and Vec2 are paired species abundance vectors for the two populations
#Adjustment is the number of clones per timepoint considered to overlap if there are 0 overlapping clones. Default 1, so 1^2 overlap. 
#Also if there are no clones in either timepoint, there is considered to be a clone group of size (1+Adjustment) at each timepoint.
MorisitaDistance <- function(Vec1,Vec2, AsProportion, Adjustment, doLogit, useHorn = FALSE) {
  Vec1 <- as.numeric(Vec1)
  Vec2 <- as.numeric(Vec2)
  Tot1 <- sum(Vec1)
  Tot2 <- sum(Vec2)
  
  SumBoth <- 0
  SumA <- 0
  SumB <- 0
  
  for (i in 1:length(Vec1)) {
    if (useHorn) {
      SumA = SumA + Vec1[i]*(Vec1[i])/(Tot1*(Tot1))
      SumB = SumB + Vec2[i]*(Vec2[i])/(Tot2*(Tot2))
    }
    else {
      SumA = SumA + Vec1[i]*(Vec1[i]-1)/(Tot1*(Tot1-1))
      SumB = SumB + Vec2[i]*(Vec2[i]-1)/(Tot2*(Tot2-1))
    }
    SumBoth = SumBoth + 2*Vec1[i]*Vec2[i]/(Tot1*Tot2)
  }
  
  if (SumBoth == 0) {
    SumBoth = 2*Adjustment^2/(Tot1*Tot2)
  }
  if (SumA+SumB == 0) {
    if (useHorn) {
      SumA <- Adjustment*(Adjustment)/(Tot1*(Tot1))
      SumB <- Adjustment*(Adjustment)/(Tot2*(Tot2))
    }
    else {
      SumA <- Adjustment*(1+Adjustment)/(Tot1*(Tot1-1))
      SumB <- Adjustment*(1+Adjustment)/(Tot2*(Tot2-1))
    }
  }
  memory <- c(SumA, SumB, SumBoth)
  ScaleFactor <- 100/(SumBoth+SumA+SumB)
  SumBoth <- SumBoth*ScaleFactor
  SumA <- SumA*ScaleFactor
  SumB <- SumB*ScaleFactor
  
  out <- 1-SumBoth/(SumA+SumB)
  
  if(AsProportion) {
    out <- SumA + SumB
  }
  
  if (doLogit == TRUE) {
    out <-  log(out/(1-out))
  }
  
  list(a=out, b=c(SumA,SumB,SumBoth), c = memory)
}


#Function that applies the Morisita distance to each pair of populations in a set of populations of clones.
#Arguments are the same as arguments to MorisitaDistance (with simplified clones table, so more than 1 row, and no rows specified).
# Output is a distance matrix
getDistanceMatrix <- function(ClonesArray, AsProportion, Adjustment, DivideBy100AndLogit, useHorn = FALSE) {
  OutputMatrix <- matrix(0, nrow(ClonesArray), nrow(ClonesArray)) 
  rownames(OutputMatrix) <- rownames(ClonesArray)
  colnames(OutputMatrix) <- rownames(ClonesArray)
  for (i in 1:nrow(OutputMatrix)) {
    for (j in 1:ncol(OutputMatrix)) {
      OutputMatrix[i,j] <- MorisitaDistance(ClonesArray[i,], ClonesArray[j,],
                                            AsProportion, Adjustment, FALSE, useHorn = useHorn)$a
      if(DivideBy100AndLogit) {
        OutputMatrix[i,j] <- max(log((OutputMatrix[i,j]/100)/(1 - OutputMatrix[i,j]/100)), 0)
      }
    }
  }
  as.dist(OutputMatrix)
}


#Vec1 and Vec2 are paired species abundance vectors for the two populations
#NumOverlap returns a list with the number of (not necessarily unique) clones that are shared at least once at each timepoint, then the number just at the first timepoint, then just at the second timepoint
NumOverlap <- function(Vec1, Vec2) {
  Vec1 <- as.numeric(Vec1)
  Vec2 <- as.numeric(Vec2)
  Tot1 <- sum(Vec1)
  Tot2 <- sum(Vec2)
  SumA <- 0
  SumB <- 0
  SumBoth <- 0
  for (i in 1:length(Vec1)) {
    if (Vec1[i]> 0 && Vec2[i] > 0) {
      SumBoth <- SumBoth + Vec1[i] + Vec2[i]
    }
    else {
      if (Vec1[i] > 1) {
        SumA <- SumA + Vec1[i]
      }
      if (Vec2[i] > 1) {
        SumB <- SumB + Vec2[i]
      }
    }
  }
  list(a = SumBoth, b = SumA, c=SumB)
}



###### HYPOTHESIS TESTS, RESAMPLING.



#There are a few hypothesis tests.
#With a single pair of timepoints, we can use a permutation test to see if Morisita center is different from 50%. THIS TESTS THE SHARP NULL that samples are from the same timepoint, but MORISITA = 50% IS RAPIDLY ASYMPTOTICALLY EQUIVALENT TO THE SHARP NULL UP TO SOME FACTOR OF o(1) DEVIATION (although I haven't proved this or what an o(1) deviation means here). 


#IMPORTANT: When we do test on TCRs, we also apply permutation to the entire process starting from replicates. Thus, we need the original info in Supp Data 4. We permute whetther each replicate belongs to timepoint 1 or timepoint 2, and then recalculate clone sizes with Poisson correction, before recalculating Morisita. 

#With paired pairs of timepoints, we can test if Morisita difference between one pair of timepoints is different from the other. There is no clear sharp null here, which makes permutation tests hard. However, we try one, along with following tests:

# 1. Delta method - z-test: Assume that variance is small so that Morisita has only small variations about its expectation. Then these variations of Morisita are approximately a linear functon of variations in clone sizes, thus the variance of Morisita is these squared linear factors times the variance of clone sizes.

#Because Morisita is a simple function of second order U-statistics it is asymptotically normal (i.e. the Simpson index and Morisita numerator are U-statistics, mildly weighted by whether sizes are from clones 1 or 2). Then, we can use a z-test between two Morisita centers, where the variance of Morisita1-Morisita2 is the sum of the Delta-estimated variances of each Morisita index, and the expectation of Morisita1-Morisita2 is of course 0 under the null.

# 2. Chi-squared method with n = number clones: Morisita essentially is sorting each pair of clones into being from the first timepoint (left flank), second timepoint (right flank), or both timepoints (center). Thus, if we know the "effective sample size", we could represent this as a 2x2 table with one column being clone pairs in the center, and the second column being clone pairs in the flanks, and the row sums being the effective sample sizes. Here each row is a pair of timepoints. We can then do a Chi-square test assuming things are asymptotically normal. 

# We can set the effective sample size for a pair of timepoints to be the total number of proviruses seen more than once (NOT the number of clone pairs, which is too high as clones are then counted repeatedly), which empirically works pretty well. This somewhat suffers from (1) assumptions on effective clone size, and (2) chi-square does a pooled variance estimate from Morisita1 = a/(a+b) and Morisita2 = c/(c+d), which isn't quite appropriate here, but benefits from (3) being well-known test and (4) working somewhat well empirically.
                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                    
# 3. Subtracted bootstrapping with replacement (m of n bootstrapping without replacement is possible too for large sample sizes but is incredibly slow because need to model variance of Morisita as function of m).: Because it's an analytic function of U-statistics (see a couple paragraphs up), Morisita sampling distribution should be recoverable asymptotically by bootstrapping with replacement.

#We INDEPENDENTLY bootstrap Morisita1 and Morisita2. For each one, we bootstrap by resampling clones1 with replacement and clones2 independently with replacement, and calculating the Morisita distance between the new clones1 and clones2. This yields resampled distributions of clones1 and clones2, allowing resampled distribution of Morisita1. We can then subtract bootstrapped distributions of Morisita1 and Morisita2 to get an estimate of resampled distributions of Morisita1-Morisita2. As suggested by Effron we can make an asymptotic 2-sided hypothesis test based on how much of this distribution overlaps with 0.


# For each testing option, we give the option to (1) calculate a p-value, and (2) for validation, simulate the form of the estimated null sampling distribution (which could be visually compared to the actual sampling distribution). For bootstrapping we estimate the resampled distribution, slightly different from null distribution. 



#### One-Morisita permutation test

#****IMPORTANT: When we do test on TCRs, we also apply permutation to the entire process starting from replicates. Thus, we need the original info in Supp Data 4. We permute whetther each replicate belongs to timepoint 1 or timepoint 2, and then recalculate clone sizes with Poisson correction, before recalculating Morisita. The process below however is simpler and what we do for proviruses.
#*See next code block for the TCR permutation code; it gets its own section because it has an extra argument compared to these functions.

# see description in text above.
# takes clones1 and clones2, paired vectors of clone sizes (each index is size of a clone).
# also takes simulatedPoints, the number of morisitas to simulate
# under assumption samples are from same timepoint, merges them and then re-partitions back into clones1 and clones2 and calculates Morisita between them to get null distribution.
# returns a list where first element is p-value, and second is simulated Morisita null distribution
# uses helper function morisitaRepeatedMixResample which just loops morisitaResampleWithElongated and is found in other file.
# if doPlot = TRUE, plots histogram of permuted morisitas.
permOneMorisitaTest <- function(clones1, clones2, simulatedPoints = 1000, doPlot = FALSE) {
  if (simulatedPoints == 0) return(1)
  morisita = getMorisita(clones1, clones2)
  sim = morisitaRepeatedMixResample(clones1, clones2, simulatedPoints)
  
  if (doPlot) {
    hist(sim, xlim=c(0,.8), breaks=20, col="#16cdd0")
    abline(v=morisita, col="red", lwd=5)
    abline(v=quantile(sim, 0.05), col="black", lwd = 5)
  }
  
  list(pval = mean(sim < morisita), sim = sim)
}

# returns 5'th percentile of permOneMorisitaTest permutation test and its p-value in a 2-element vector
# arguments are as in permOneMorisitaTest; no plotting done.
morisitaResamplingSignificance <- function(clones1, clones2, simulatedPoints) {
  all = permOneMorisitaTest(clones1, clones2, simulatedPoints)
  morisita = getMorisita(clones1, clones2)
  return(c(as.numeric(100*quantile(x=all$sim, probs=c(.05))), all$pval))
}

### DELTA TEST
# Uses functions in helper code file to get estimated variance with delta method.
# Now that we can estimate the variance of Morisita1-Morisita2 as demonstrated in next block of code, we can calculate a p-value and sample from the null distribution.
# clonesA1 and clonesA2 are used for morisita1, and clonesB1 and clonesB2 are used for morisita2. simulatedPoints is number of simulated draws from estimated null distribution.
# returns a p-value, and if sumulatedPoints isn't 0, a list with first element the p-value and second a numeric vector of samples from null distribution with number of samples equal to simulatedPoints.
deltaTest <- function(clonesA1, clonesA2, clonesB1, clonesB2, simulatedPoints = 0) {
  observed = getMorisita(clonesA1, clonesA2) - getMorisita(clonesB1, clonesB2)
  estimateVariance = getMorisitaDerivatives(clonesA1, clonesA2) + 
    getMorisitaDerivatives(clonesB1, clonesB2)
  
  pval = 2*pnorm(q=abs(observed/sqrt(estimateVariance)), lower.tail=FALSE)
  if (simulatedPoints == 0) return(pval)
  
  sim = rnorm(simulatedPoints, mean = 0, sd = sqrt(estimateVariance))
  
  list(pval = pval, sim = sim)
}

# Takes either simple tables or objects object1 and object2 (gets simple tables from simple objects if isObject = T)
# i1 and j1 are indices of timepoints in object1 between which you're getting a morisita1
# i2 and j2 are indices of timepoints in object1 between which you're getting a morisita2
# Returns 2-tailed z-test p-value that morisita1 and morisita2 are not equal, with variance estimated by delta method.
getDeltaZTestFromObjects <- function(object1, object2, i1, j1, i2, j2, isObject=T) {
  if (isObject) {
    object1 = object1$simpleTable
    object2 = object2$simpleTable
  }
  morDiffEst <- getMorisita(object1[i1,], object1[j1,])-getMorisita(object2[i2,], object2[j2,])
  varDiffEst <- deltaMorisitaVariance(as.numeric(object1[i1,]), as.numeric(object1[j1,])) + 
    deltaMorisitaVariance(as.numeric(object2[i2,]), as.numeric(object2[j2,]))
  morDiffEst/sqrt(varDiffEst)
  2*pnorm(q=abs(morDiffEst/sqrt(varDiffEst)), lower.tail=FALSE)
}


# takes a simple object and indices of two timepoints, and isObject which if false lets you make object1 a simple table instead of simple object.
# returns a 95% confidence interval on Morisita overlap between timepoints by Delta method.
getDelta95Interval <- function(object1, i, j, isObject = T) {
  if (isObject) object1 <- object1$simpleTable
  deltaEstimatedSD1 <- sqrt(deltaMorisitaVariance(object1[i,], object1[j,]))
  return(c(max(0,100*(getMorisita(object1[i,], object1[j,]) - 1.96*deltaEstimatedSD1)),100*(getMorisita(object1[i,], object1[j,]) + 1.96*deltaEstimatedSD1)))
}


#### Chi-square test
# see description above. Makes 2x2 table where row sums are number of clones in the two timepoints used to calculate Morisita1 (row 1) and Morisita2 (row 2), 
# first column is normalized number of pairs of clones that would be sorted to center, and second is number of pairs of clones that would be sorted to flanks.
# then do chi-squared
# clonesA1 and clonesA2 are used for morisita1, and clonesB1 and clonesB2 are used for morisita2. simulatedPoints is number of simulated draws from estimated null distribution.
# returns a p-value, and if sumulatedPoints isn't 0, a list with first element the p-value and second a numeric vector of samples from null distribution with number of samples equal to simulatedPoints.
chiTest <- function(clonesA1, clonesA2, clonesB1, clonesB2, simulatedPoints = 0) {
  morisitaA = getMorisita(clonesA1, clonesA2)
  morisitaB = getMorisita(clonesB1, clonesB2)
  
  mergedA = clonesA1 + clonesA2
  nEffectiveA = sum(mergedA[mergedA>1]) # number of clones is effective sample size = row sum
  mergedB = clonesB1 + clonesB2
  nEffectiveB = sum(merged[merged>1]) 
  # 2x2 table as described:
  x = matrix(c(morisitaA*nEffectiveA, (1-morisitaA)*nEffectiveA, 
               morisitaB*nEffectiveB, (1-morisitaB)*nEffectiveB), nrow=2,ncol=2) 
  
  pval = chisq.test(x)$p.value
  if (simulatedPoints == 0) return(pval)
  
  presumedProb = (morisitaA*nEffectiveA + morisitaB*nEffectiveB)/(nEffectiveA + nEffectiveB) # pooled Morisita estimated
  estimateVariance = presumedProb*(1-presumedProb)*(1/nEffectiveA + 1/nEffectiveA) # pooled variance estimate
  sim = rnorm(simulatedPoints, mean = 0, sd = sqrt(estimateVariance))
  
  
  list(pval = pval, sim = sim)
}


# Takes either simple tables or objects object1 and object2 (gets simple tables from simple objects if isObject = T)
# i1 and j1 are indices of timepoints in object1 between which you're getting a morisita1
# i2 and j2 are indices of timepoints in object1 between which you're getting a morisita2
# Returns chi-square p-value that morisita1 and morisita2 are not equal, with sample size estimated as number of repeated sequences.
# Note this is equivalent to a pooled-variance two-proportion Z-test where variance is calculated as if number of repeated objects is number of repeated sequences
chiSquareTestFromObjects <- function(object1, object2, i1, j1, i2, j2, isObject = T) {
  if (isObject) {
    object1 = object1$simpleTable
    object2 = object2$simpleTable
  }
  morisita1 = getMorisita(object1[i1,], object1[j1,])
  morisita2 = getMorisita(object2[i2,], object2[j2,])
  x=object1[i1,] + object1[j1,]
  n1 = sum(x[x>1])
  x=object2[i2,] + object2[j2,]
  n2 = sum(x[x>1]) 
  x = matrix(c(morisita1*n1, (1-morisita1)*n1, 
               morisita2*n2, (1-morisita2)*n2), nrow=2,ncol=2)
  round(chisq.test(x)$p.value,8)
}


#### Double permutation test
# see description at top.
# Calculates permutation distribution of morisita for each pair of timepoints like in permOneMorisitaTest, but then subtracts them to get a null distribution of difference in Morisitas.
# clonesA1 and clonesA2 are used for morisita1, and clonesB1 and clonesB2 are used for morisita2. simulatedPoints is number of simulated draws from estimated permutation distribution.
# returns a list with first element the p-value and second a numeric vector of samples from null distribution with number of samples equal to simulatedPoints.
permTwoMorisitaTest <- function(clonesA1, clonesA2, clonesB1, clonesB2, simulatedPoints = 1000) {
  if (simulatedPoints == 0) return(1)
  observed = getMorisita(clonesA1, clonesA2) - getMorisita(clonesB1, clonesB2)
  
  sim = morisitaRepeatedMixResample(clonesA1, clonesA2, simulatedPoints) - 
    morisitaRepeatedMixResample(clonesB1, clonesB2, simulatedPoints)
  
  # two-sided pval as double min of one-sided pvals:
  pval = 2*min(mean(sim <= observed), mean(sim >= observed)) 
  
  list(pval = pval, sim = sim)
}


#### Bootstrap test
# relies on helper getBootstrapDistribution function in helper function file to bootstrap, which is very simple: resample clone sizes at all timepoints with replacement, and calculate Morisita between resampled clone sizes.
# clonesA1 and clonesA2 are used for morisita1, and clonesB1 and clonesB2 are used for morisita2. simulatedPoints is number of simulated draws from estimated bootstrap distribution.
# returns a list with first element the p-value and second a numeric vector of samples from bootstrap distribution with number of samples equal to simulatedPoints. 
# NOTE THIS IS ONLY CASE WHERE RETURNED DISTRIBUTION IS NOT A NULL DISTRIBUTION.
bootstrapTest <- function(clonesA1, clonesA2, clonesB1, clonesB2, simulatedPoints = 1000) {
  
  if (simulatedPoints == 0) return(1)
  
  morisitaA = getMorisita(clonesA1, clonesA2)
  morisitaB = getMorisita(clonesB1, clonesB2)
  observed = morisitaA - morisitaB
  
  
  sim = getBootstrapDistribution(clonesA1, clonesA2, nBoots = simulatedPoints) - 
    getBootstrapDistribution(clonesB1, clonesB2, nBoots = simulatedPoints)
  sim = sim + (observed - mean(sim)) # center the bootstrap on observed value.
  
  # two-sided pval as double min of one-sided pvals:
  pval = 2*min(mean(sim <= 0), mean(sim >= 0)) 
  
  list(pval = pval, sim = sim)
}


# Samples permutation distribution for TCR morisita
# Takes the full matrix of supp data 4 containing each column as a replicate each row as a clone
# and each cell as the number of RNA reflecting the clone  at the replicate.
# Also takes the year labels of two years (in Supp Data 4 these are years on ART whereas for the RDS above they're absolute years, but just need to be consistent to run the function).
# There's only a limited number of possible unique permutations (e.g. if 12 replicates at two timepoints, there's 24 choose 12 possible label assignments); these are sampled without replacement, and if maxResamples < number possible permutations it uses all of them. If there's more than 24 total replicates we just sample permutations with replacement for simplicity since enumerating them all would take too much space (could go back and eliminate duplicate permutations in theory)
# Returns a list where first element is a numeric vector of length 1 with the true Morisita and the second element is a numeric vector of all the permutation-sampled morisitas.

TCRResamplingsTwoTimepoints <- function(mat, year1, year2, maxResamples=2000) {
  
  ### some matrix cleanup since we just want the "Rep" columns and we don't want rows for Simpson and PerClones which don't reflect actual clone sizes. Also sometimes we did 12 replicates twice at a timepoint and the second set are marked "B"; we want to merge them.
  mat <- mat[,grep("Rep", colnames(mat))]
  mat <- mat[,-grep("NumPositiveReps", colnames(mat))]
  
  if (length(grep("Simpsons", rownames(mat))) > 0) mat <- mat[-grep("Simpsons", rownames(mat)),]
  if (length(grep("PercClones", rownames(mat))) > 0) mat <- mat[-grep("PercClones", rownames(mat)),]
  
  
  
  Cols1 <- which(sub("\\-.*", "", colnames(mat)) %in% c(year1, paste0(year1, "B")))
  Cols2 <- which(sub("\\-.*", "", colnames(mat)) %in% c(year2, paste0(year2, "B")))
  
  # Only focus on the relevant cols.
  mat2 <- data.frame(mat[,c(Cols1, Cols2)])
  
  # Use presence-absence at each timepoint rather than RNA abundance.
  for (j in 1:ncol(mat2)) {
    mat2[,j] <- as.numeric(mat2[,j])
    mat2[,j] <- ifelse(mat2[,j]>0,1,0)
  }
  
  
  # We make permuts, which is a matrix where each row is for one permutation and contains the indices of the columns of mat2 which are assigned to timepoint 1 in the permutation.
  # if less or equal to 24 choose 12 permutations, we use combn to make all the possible permutations of columns that could be assigned to timepoint 1 
  # For example, if replicates 1-12 are time 1 and replicates 13-24 are time 2, a row of permut might be 1,4,5,6,14,15,16,18,19,20,21,24, representing replicates assigned to time 1 in a permutation
  # The rows of this matrix are then subsampled without replacement to get maxResamples final permutations.
  # If more than 24 choose 12 possible permutations, enumerating everything would make too big a matrix, so we use gtools::permute function to randomly sample a permutation for each row. In theory rows could be identical, but should be quite rare - I believe could theoretically check for duplicate rows in O(maxResamples) time with a hash map.
  if (ncol(mat2) <= 24) {
    permuts=t(combn(1:ncol(mat2), length(Cols1)))
    if (nrow(permuts) > maxResamples) {
      permuts <- permuts[sample.int(nrow(permuts), maxResamples),]
    }
  }
  
  else {
    permuts <- matrix(0, maxResamples, length(Cols1)) 
    for (i in 1:maxResamples) {
      permuts[i,] <- gtools::permute(1:ncol(mat2))[1:length(Cols1)]
    }
  }
  
  # Other is simply the indices of replicates which are assigned to timepoint 2 in each permutation; they're just the leftover indices from each row of permut.
  other <- t(apply(permuts, 1, function(x){c(1:ncol(mat2))[-which(c(1:ncol(mat2))%in% x)]}))
  
  
  # Now for each replicate permutation, we sum the number of positive replicates and Poisson-correct it for the permutation timepoint1 and timepoint2 replicates (derived from a matched row of permuts and other matrices). This gives us clone sizes for permuted timepoints which can be compared.
  # This is done by the morisitaFromReadCols function in the helper functio file.
  mat2 <- as.matrix(mat2)
  mat2 <- data.frame(mat2)
  outs <- numeric(nrow(permuts))
  for (i in c(1:nrow(permuts))) {
    outs[i] <- morisitaFromReadCols(mat2, permuts[i,], other[i,])
  }
  
  # also gets the original Morisita, to return as list.
  original = morisitaFromReadCols(mat2, 1:length(Cols1), (length(Cols1)+1):ncol(mat2))
  
  list(original = original, sim = outs)
  
}





############### Codes that help with hypothesis tests

#### DELTA METHOD
# takes clones1 and clones2, paired vectors of clone sizes (each index is size of a clone).
# calculates the derivatives of Morisita center at observed clone sizes, and uses this to estimate variance of Morisita.
# similar to mueller and altenberg 1985 except derivatives slightly alterred since we use morisita center instead of flanks divided by center.
deltaMorisitaVariance <- function(clones1, clones2) {
  clones1 = as.numeric(clones1)
  clones2 = as.numeric(clones2)
  derivativesOut <- getMorisitaDerivatives(clones1, clones2)
  estimatedProportions <- estimateMultinomialProportions(clones1, clones2, derivativesOut$goodIndex) #esssentially clones1/sum(clones2) and clones2/sum(clones2)
  # Now we just find x^TVx where x is the gradient and V is the multinomial covariance, to estimate variance of Morisita.
  # However, since V is rank 1, the calculation below is equivalent and avoids storing in memory the massive covariance matrix.
  # Also note variance due to clones1 variance is independent from that due to clones2 variance, so variance contributions of each can be calculated separately and summed.
  varianceEst1 <- ((derivativesOut$der1)^2 %*% estimatedProportions[[1]] - (derivativesOut$der1%*%estimatedProportions[[1]])^2)/sum(clones1)
  varianceEst2 <- ((derivativesOut$der2)^2 %*% estimatedProportions[[2]] - (derivativesOut$der2%*%estimatedProportions[[2]])^2)/sum(clones2)
  as.numeric(varianceEst1+varianceEst2)  # adds variance due to independent variances of clones1 and clones2 to get variance of Morisita.
}

#more Delta helper functions:

# gets derivative vector (gradient) of morisita center with respect to each clone size (or technically the probability of each clone, which is clones1/sum(clones1)), calculated at the observed clone sizes (which is a plug-in estimate). 
# technically the last clone size is not random at fixed sample size, so to avoid a singular covariance matrix we consider it to be fixed (which is just the clone at goodIndex)
# similar to mueller and altenberg 1985 except derivatives slightly alterred since we use morisita center instead of flanks divided by center.
# returns a list where the first element is the index of the fixed-size clone, the second element is the gradient of Morisita at the observed clone sizes of clones1, and the third element is the gradient of Morisita at the observed clone sizes of clones2.
getMorisitaDerivatives <- function(clones1, clones2) {
  prop1 <- as.numeric(clones1/sum(clones1))
  prop2 <- as.numeric(clones2/sum(clones2))
  
  goodIndex <- which(prop1 > 0 & prop2 > 0)[1]
  
  pre1 <- prop1[-goodIndex]
  pre2 <- prop2[-goodIndex]
  morisitaDenominator <- getMorisitaDenominator(clones1, clones2)
  morisitaNumerator <- getMorisitaNumerator(clones1, clones2)
  der1 <- 2*(pre2-prop2[goodIndex] - (morisitaNumerator/morisitaDenominator)*(pre1-prop1[goodIndex] + pre2-prop2[goodIndex]))/morisitaDenominator
  der2 <- 2*(pre1-prop1[goodIndex] - (morisitaNumerator/morisitaDenominator)*(pre2-prop2[goodIndex] + pre1-prop1[goodIndex]))/morisitaDenominator
  
  list(goodIndex=goodIndex, der1=der1, der2=der2)
}

# simple function that just gives plug-in estimates of multinomial proportions of all clones except for the index of clone we consider to have fixed size.
estimateMultinomialProportions <- function(clones1, clones2, goodIndex) {
  sum1 <- sum(clones1)
  sum2 <- sum(clones2)
  prop1t <- clones1[-goodIndex]/sum1
  prop2t <- clones2[-goodIndex]/sum2
  list(prop1t, prop2t)
}

#### PERMUTATION
### Does the permutation morisitas for permutation tests.
morisitaRepeatedMixResample <- function(clones1, clones2, nResamples) {
  merged <- clones1+clones2
  elongatedForm <- elongate(merged)
  all <- numeric(nResamples)
  for (i in 1:nResamples) {
    all[i] <- morisitaResampleWithElongated(clones1, clones2, merged, elongatedForm)
  }
  all
}

### BOOTSTRAPPING
# takes numeric vectors of clone sizes clones1 and clones2, and number of bootstraps to do nBoots.
# resamples clones1 and clones2 with replacement, i.e. assuming multinomial proportions equal to plug-in estimate of clones1/sum(clones1) and clones2/sum(clones2).
# then calculates distribution of morisitas between resampled clones1 and clones2.
getBootstrapDistribution <- function(clones1, clones2, nBoots) {
  real1 <- clones1/sum(clones1)
  real2 <- clones2/sum(clones2)
  
  resamples1 <- rmultinom(nBoots, sum(clones1), real1)
  resamples2 <- rmultinom(nBoots, sum(clones2), real2)
  bootMors <- numeric(nBoots)
  for (i in 1:nBoots) {
    bootMors[i] <- getMorisita(resamples1[,i], resamples2[,i])
  }
  bootMors
}


# helps with subsampling.
generateSubsample = function(vec, sampSize, replace = FALSE, elongated = NULL) {
  out = list()
  len = length(vec)
  if (sum(vec) < 2E8 && replace == FALSE) {
    if (is.null(elongated)) elongated = elongate(vec)
    resample <- tabulate(sample(elongated, sampSize, FALSE))
    if (len > length(resample)) resample <- c(resample, rep(0, len-length(resample)))
  }
  else {
    resample <- tabulate(sample(1:length(vec), sampSize, TRUE, prob = vec))
    if (len > length(resample)) resample <- c(resample, rep(0, len-length(resample)))
  }
  resample
}

# helps with quantiles.
getPQuantiles <- function(dist, p) {
  c(quantile(dist, p/2),quantile(dist, 1-p/2))
}




######## PERMANOVA helping funcitons

logitFunc = function(x) log(x/(1-x))

#makes "distance"/similarity matrix of Morisita similarity between each pair of rows in a simple table of a simple Object. Can make them logits of the similarity if we want.
getMorisitaSimilarities <- function(simpleCloneObject, useLogit = TRUE) {
  simpleCloneTable = simpleCloneObject$simpleTable
  simpleCloneObject$newTimes = as.numeric(simpleCloneObject$newTimes)
  outSims <- outTimes <- matrix(0, nrow(simpleCloneTable), nrow(simpleCloneTable))
  for (i in 1:nrow(outSims)) {
    for (j in 1:ncol(outSims)) {
      morisita = getMorisita(simpleCloneTable[i,], simpleCloneTable[j,], avoid0=T)
      outSims[i,j] = ifelse(useLogit, logitFunc(morisita), morisita)
      outTimes[i,j] = abs(simpleCloneObject$newTimes[i] - simpleCloneObject$newTimes[j])
    }
  }
  list(sims = as.dist(outSims), times = as.dist(outTimes))
}


######## OTHER

# takes a species abundance vectors for a single population, and returns the proportion of organisms belonging to a species present more than once.
percentClones <- function(vec) {
  vec = as.numeric(vec)
  sum(vec[vec>1])/sum(vec)
}
