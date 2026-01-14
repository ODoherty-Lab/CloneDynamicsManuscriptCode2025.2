# CloneDynamicsManuscriptCode2025
Code for manuscript in review on HIV-infected T cell clonal dynamics. 

This will be updated with more detailed instructions. In particular, supplementary data numbers will change. However, hopefully should be able to tell from the titles of supplementary data what the appropriate names are. Please reach out if any element of the code is unclear or not running as might be expected. As run locally, this code produces the figures and tables in the paper from raw data. 

Note that for the Morisita data, can start by running SuppDataToRDSObjects.R on Supp Data to reproduce RDS objects as in RDS objects folder. These are then used in ProduceMorisitaOutputs.Rmd which walks through each figure. 


Functions are defined in the .R files of the HelpingFunctions folder. Simple objects in the Morisita code parts are just 2-element lists with first element is "simpleTable" and second is "newTimes".
simpleTable: a numeric data frame where a column is a clone (colnames are just numbered indices) and each row is a timepoint, gives clone size at time.
newTimes: a vector of the years, corresponding to rows of simpleTable.

Supp Data containing a clean version of source data should be downloaded from paper and added to the Supp Data folder (particular the TCR supp data was too large to add to the repository, although it is referenced in the code). Note the .RProj isn't strictly necessary but it avoids having to change any directory paths as its directory is automatically accessed by here::here() in R (should work even if you don't access the code from the .RProj file; just having it in directory is in theory enough). You may need to change the base path variable as commented at the top of the R Markdown file you work with, to help redirect code to the proper folder.

Depending on demand and constraints we might add more general information about code objects to this README page. That said, we intend the present comments, currently present throughout the code, to be  sufficient to follow all objects and functions throughout the code, so please reach out if they fail to do. 

Please be aware that supplemental data changes underwent superficial numbering changes recently, which are now reflected in the code, which should not interfere with functioning but this will be confirmed through re-testing this weekend (on 1/17/26).
