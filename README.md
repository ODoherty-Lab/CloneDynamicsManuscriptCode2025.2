# CloneDynamicsManuscriptCode2025
Code for manuscript in review on HIV-infected T cell clonal dynamics

This will be updated with more detailed instructions. 

Repository currently contains code to reproduce figures 2-5 and Supp Data 7-8, from Supp Data 4 and 5. Also contains helpful functions for hypotheses testing, testing hypothesis tests, and working with the data.
Labeling may differ from final figures due to e.g. header changes/column collapse in Excel, but figure portrayal of data is exactly as used in paper.

First, SuppDataToRDSObjects.R can be run on Supp Data to reproduce RDS objects as in RDS objects folder. These are then used in ProduceMorisitaOutputs.Rmd which walks through each figure. 

Functions are defined in the .R files of the HelpingFunctions folder. Simple objects are just 2-element lists with first element is "simpleTable" and second is "newTimes".
simpleTable: a numeric data frame where a column is a clone (colnames are just numbered indices) and each row is a timepoint, gives clone size at time.
newTimes: a vector of the years, corresponding to rows of simpleTable.

Supp Data 4 and 5, containing a clean version of source data, should be downloaded from paper and added to the Supp Data folder.
