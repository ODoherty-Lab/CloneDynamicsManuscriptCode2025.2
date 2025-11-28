

# takes a "tabulatedVersion" list. 
# This list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.


# cutoffs gives the clone size cutoffs on the columns of the proportion sense plots.
# cutoffs is a numeric vector of any length > 1. Each value is the inclusive cutoff of one of the bars on the proportion plot.
# e.g. if cutoffs are c(1,2,4,7), you'll get bars for clone size >=1, size >=2, >=4, and >= 7.
# oneSided can be 0 (FALSE), 1, or 2, reflecting 2-sided, greater alternative, and lesser alternative respectively, for performing proportion plot binomial tests.

# returns a list of vectors, all the same length. The index of an element in each vector is reflective of its cutoff.
# Note this is only done for integrations in significant genes (e.g. in growth-related genes), reflected in the first
# two vectors of "tabulatedVersion" which are from the rows in the original inTable where the 4th column is TRUE.
# Vectors in the list provide the proportion sense, the binomial p-value of the %sense being not equal to 50%, 
# the lower CIs and upper CIs on the proportion sense, and the cumulative inclusive clone size cutoffs associated with each index of the list's vectors.

getProportionSenseAtCutoffsWithBinomials <- function(tabulatedVersion, cutoffs = NULL, oneSided = FALSE) {
  
  tabulatedVersion <- padVecs(tabulatedVersion) # Makes all four lists the same length by adding 0's
  # specifically, this is saying that even if no sense integrations in significant genes of size 40 might be observed,
  # I make my vector of clone sizes stretch out to 40, noting that there's 0 sense significant clones of size 40.
  
  ObsCPlus <- tabulatedVersion[[1]]  # clone sizes of sense integrations into significant (e.g. growth-related) genes. The i'th element is the number of distinct size-i clones.
  ObsCMinus <- tabulatedVersion[[2]] # clone sizes of antisense integrations into significant (e.g. growth-related) genes. The i'th element is the number of distinct size-i clones.
  fullLength <- length(ObsCPlus) # note that after padding, this is just the size of the largest observed clone.
  
  # This produces a vector of clone sizes of all integrations in significant genes 
  # note untabulate transforms from a tabulatedVersion vector (where the i'th element is the number of clones of size i)
  # into a vector of clone sizes, where in the untabulated version, if "3" occurs 4 times, then there's 4 distinct clones of size 3.
  values <- c(untabulate(ObsCPlus), untabulate(ObsCMinus)) 
  
  
  # catch empty tests.
  if (length(values) == 0) {
    return(list(proportions = c(), binomials=c(), cutoffsUsed = c()))
  }
  if (length(cutoffs) == 0 || fullLength == 0 || sum(ObsCPlus + ObsCMinus) <= 0) {
    return(list(proportions = c(), binomials=c(), cutoffsUsed = c()))
  }
  
  
  # remove cutoffs that would lead to proportion sense being 0/0
  highestPossibleCutoff <- max(values) # largest observed clone size.
  cutoffsToDiscard <- which(cutoffs > highestPossibleCutoff)
  cutoffsUsed <- cutoffs
  if (length(cutoffsToDiscard) > 0) cutoffsUsed <- cutoffs[-cutoffsToDiscard]
  
  if (length(cutoffsUsed) == 0) {
    return(list(proportions = c(), binomials=c(), cutoffsUsed = c()))
  }
  
  
  # prepare to output proportion of clones in sense above each cutoff, its associated p-value, and the confidence interval.
  proportions <- binomials <- lowerCIs <- upperCIs <- Ns <- numeric(length(cutoffsUsed))
  for (i in 1:length(cutoffsUsed)) { # for each cumulative cutoff
    senseNum <- sum(ObsCPlus[cutoffsUsed[i]:fullLength]) # Number of distinct clones with sense HIV in significant (e.g. growth-related) genes with size clone size at least equal to the cutoff
    antiNum <- sum(ObsCMinus[cutoffsUsed[i]:fullLength]) # Number of distinct clones with antisense HIV in significant (e.g. growth-related) genes with size clone size at least equal to the cutoff
    proportions[i] <- senseNum/(senseNum + antiNum) # proportion sense among significant clones meeting the cumulative inclusive cutoff.
    tester <- binom.test(senseNum, senseNum + antiNum) # binomial test of this
    binomials[i] <- tester$p.value # associated p-value
    lowerCIs[i] <- tester$conf.int[1] # associated lower CI
    upperCIs[i] <- tester$conf.int[2]# associated upper CI
    if (oneSided == 1) { # in this case redo the test with the "greater" alternative
      tester <- binom.test(senseNum, senseNum + antiNum, alternative = "greater")
      binomials[i] <- tester$p.value
      lowerCIs[i] <- tester$conf.int[1]
      upperCIs[i] <- tester$conf.int[2]
    }
    if (oneSided == 2) {# in this case redo the test with the "lesser" alternative
      tester <- binom.test(senseNum, senseNum + antiNum, alternative = "less")
      binomials[i] <- tester$p.value
      lowerCIs[i] <- tester$conf.int[1]
      upperCIs[i] <- tester$conf.int[2]
    }
    
  }
  
  list(proportions = proportions, binomials=binomials, cutoffsUsed = cutoffsUsed,
       lowerCIs = lowerCIs, upperCIs = upperCIs)
}



# takes a list of numeric vectors.
# Makes them all the same length by appending 0s to all but the smallest vector in the list.
# Returns the list of these padded vectors.
padVecs <- function(vecs) {
  totalLength <- max(sapply(vecs, length))
  for (i in 1:length(vecs)) {
    if (length(vecs[[i]]) < totalLength) {
      vecs[[i]] <- c(vecs[[i]], rep(0, totalLength - length(vecs[[i]])))
    }
  }
  vecs
}

# takes a numeric "tabulated" vector where the i'th element is the number of distinct clones of size i.
# Returns a vector of clone sizes which would tabulate to this tabulated vector
# For example, if tabulatedVec[4] = 3, then there's 3 size-5 clones. 
# In this case, I would add three 5's to the output vector, reflecting that there's 3 clones of size 5.
# This interpolates between the inTable sizes and the tabulatedVersion vectors.

untabulate <- function(tabulatedVec) {
  outVec <- c()
  for (i in 1:length(tabulatedVec)) {
    if (tabulatedVec[i] > 0) {
      outVec <- c(outVec, rep(i, tabulatedVec[i]))
    }
  }
  outVec
}



# cumulative logit regression test.
# Note this is restricted integrations in significant genes (e.g. in growth-related genes), reflected in the first
# two vectors of "tabulatedVersion" which are from the rows in the original inTable where the 4th column is TRUE.

# Recall tabulatedVersion is a list of four elements. Each of these is a vector where the i'th index reflects the number of distinct clones of size i.
# The first element of the list reflects sense clones in significant genes (i.e. highly-expressed or growth-related, as related in CancerRow column of table).
# The second element reflects antisense clones in significant genes
# The third and fourth elements reflect sense and antisense clones, respectively, in nonsignificant genes.

# regresses the clone size on the orientation using the polr function in R. Returns vector of 3 elements, the first two being regression coefficients and last being p-value.
performPolrTest <- function(tabulatedVersion) {
  
  # note untabulate transforms from a tabulatedVersion vector (where the i'th element is the number of clones of size i)
  # into a vector of clone sizes, where in the untabulated version, if "3" occurs 4 times, then there's 4 distinct clones of size 3.
  # then length(untabulate(x)) = sum(x) = number of distinct clones.
  # The resulting data frame here therefore has each row corresponding to a clone, and the first column is its size, and second is its orientation.
  activeOnlyContrast <- as.data.frame(cbind(c(untabulate(tabulatedVersion[[1]]), untabulate(tabulatedVersion[[2]])), 
                                            c(rep(1, sum(tabulatedVersion[[1]])), rep(0, sum(tabulatedVersion[[2]])))))
  colnames(activeOnlyContrast) <- c("Size", "isSense") 
  #first column is size of clone in active gene, second column is 1 if sense, 0 if antisense
  # 
  
  activeOnlyContrast[,1] <- as.factor(activeOnlyContrast[,1])
  if (length(unique(activeOnlyContrast[,1])) < 3) { # can't do if not enough values.
    return(NULL)
  }
  
  a <- polr(formula = Size ~ isSense, data = activeOnlyContrast, Hess = TRUE) # runs regression.
  ctable <- coef(summary(a))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2 # p-value
  c(ctable[1,1:2], p[1])
}