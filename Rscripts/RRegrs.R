# ======================================================================
# RRegrs - R Regression Models
#
# Get the best regression models for one dataset using R caret methods
# Developed as tool for nano-toxicity QSAR models
# NTUA and UM groups, enanomapper.net
# ----------------------------------------------------------------------
# Contact: 
# Cristian R Munteanu, BiGCaT - UM, muntisa [at] gmail [dot] com
# Georgia Tsiliki, ChemEng - NTUA, g_tsiliki [at] hotmail [dot] com
# ======================================================================

# Main input file: CSV
# ----------------------------------------------------------------------
# Variable names
# -----------------------------------------------------------------------
# DataSet = dataset including header, Y and X features values
# Yname, Xnames = names of dependent and features
# NoX = no of features
# NoCases = no of cases
# YData = values of Y
# XData = values of the features
# -----------------------------------------------------------------------
# Note:
# - Input file = CSV format and extension:
#    - there is only one input file with the original dataset:
#      Dependent Variable, Feature 1, ..., Feature n
#    - the names of the variables are located in the first row
# - the main script is importing separate file functions
# - the functions are reading datafiles in specific folders
# and it is writing the output as variables or datafiles
# - it generates several dataset files for each step
# - it generates PNG plots for the correlation matrix
# before and after the correlation removal
# - all the input and output files are placed into the same folder

#==========================================================================================
# (1) Load dataset and parameters
#     (these parameters will be read from an input file! TO BE IMPLEMENTED at the end)
#==========================================================================================
# (1.1) PARAMETERS
# -----------------------------------------------------------------------
# Option to run any step
# -----------------------------------------------------------------------
fDet         = TRUE  # flag to calculate and print details for all the functions
fFilters     = TRUE  # flag to apply filters                          (2)
fScaling     = TRUE  # flag for dataset Scaling                       (3)
fRemNear0Var = TRUE  # flag for Removal of near zero variance columns (4)
fRemCorr     = TRUE  # flag for Removal of correlated columns         (5)
fFeatureSel  = FALSE  # flag for wrapper methods for feature selection (7)

cutoff       = 0.9   # cut off for correlated features
fGLM         = TRUE  # flag to run GLM (8.2)

# ----------------------------------------------------------------------------------------
iScaling = 1 # 1 = normalization; 2 = standardization, 3 = other; any other: no scaling
iScalCol = 1 # 1 = including dependent variable in scaling; 2: only all features; etc.
# ----------------------------------------------------------------------------------------
trainFrac  = 3/4 # the fraction of training set from the entire dataset
                 # 1 - trainFrac = the rest of dataset, the test set
iSplitTimes = 3 # default is 10; time to split the data in train and test (steps 6-11); report each step + average

# -------------------------------------------------------------------------------------------------------
# Files
# -------------------------------------------------------------------------------------------------------
PathDataSet    = "DataResults"            # dataset folder for input and output files
DataFileName   = "ds.csv"                 # input step 1 = ds original file name
No0NearVarFile = "ds3.No0Var.csv"         # output step 3 = ds without zero near vars
ScaledFile     = "ds4.scaled.csv"         # output step 4 = scaled ds file name (in the same folder)
NoCorrFile     = "ds5.scaled.NoCorrs.csv" # output step 5 = dataset after correction removal
ResFile        = "RRegresRezults.txt"     # the main output file

glmFile        = "8.2.GLM.details.txt"

# Generate path + file name = original dataset
inFile <- file.path(PathDataSet, DataFileName)
# Main result file (append data using: sink(outRegrFile, append = TRUE) !!!)
outRegrFile <- file.path(PathDataSet,ResFile) # the same folder as the input 


cat("======================================================================
RRegrs - R Regression Models
Get the best regression models for one dataset using R caret methods
NTUA and UM groups, enanomapper.net

Contacts:
Cristian R Munteanu - muntisa [at] gmail [dot] com
Georgia Tsiliki - g_tsiliki [at] hotmail [dot] com
======================================================================\n")

# -----------------------------------
# (1.2) Load the ORIGINAL DATASET
# -----------------------------------
cat("-> [1] Loading original dataset ...\n")
# (it can contain errors, correlations, near zero variance columns)
ds.dat0 <- read.csv(inFile,header=T)                              # original dataset frame

# resolving the text to number errors for future calculations
ds.indx<- colnames(ds.dat0)[2:dim(ds.dat0)[2]]                    # FEATURE names (no dependent variable)
ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],2:dim(ds.dat0)[2]]            # dataset as columns
ds.dat1<- apply(ds.dat1,1,function(x)as.numeric(as.character(x))) # dataset as row vectors to be used with caret!!!

# dependent variable
net.c<- ds.dat0[,1]
net.c<- as.numeric(as.character(net.c)) # values
# full ds frame with training and test
ds<- as.data.frame(cbind(net.c,t(ds.dat1)))

#========================================================
# (2) FILTERS
#     (it will be implemented in the future versions)
#========================================================
# 2.1 Outlier removal
# 2.2 Custom filter (percentage threshold)
# 2.3 Processing of missing values - use of preProcess();
#     caret employs knnImpute algorithm to impute values from a neighborhood of k
cat("-> [2] Filtering dataset ... No filter!\n")

# -----------------------------------------------------------------------
# (3) Remove near zero variance columns
# -----------------------------------------------------------------------
if (fRemNear0Var==TRUE) {
  cat("-> [3] Removal of near zero variance columns ...\n")
  outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
  
  # get the ds without near zero cols 
  source("s3.RemNearZeroVar.R")                    # add function         
  ds <- cbind("net.c" = ds[,1],RemNear0VarCols(ds[,2:dim(ds)[2]],fDet,outFile))
  # use df without Y (predicted values), reconstruct the ds
  # inputs: ds, flag for details, output file
}

# -----------------------------------------------------------------------
# (4) Scaling dataset: normalization (default), standardization, other
# -----------------------------------------------------------------------
if (fScaling==TRUE) {
  cat("-> [4] Scaling original dataset ...\n")
  outFile <- file.path(PathDataSet,ScaledFile)       # the same folder as input
  
  # run fuction for scaling input dataset file
  source("s4.ScalingDataSet.R")                      # add function
  ds <- ScalingDS(ds,iScaling,iScalCol,fDet,outFile)
  # use df without Y (predicted values), reconstruct the ds
  # inputs: ds, type of scaling, flag for details, starting column, output file
}

# -----------------------------------------------------------------------
# (5) Remove correlated features
# -----------------------------------------------------------------------
if (fRemCorr==TRUE) {    
  cat("-> [5] Removing correlated features ...\n") 
  outFile <- file.path(PathDataSet,NoCorrFile)    # the same folder as the input
  
  # run function to remove the correlations between the features
  source("s5.RemCorrFeats.R")                     # add function
  ds <- cbind("net.c" = ds[,1],RemCorrs(ds[,2:dim(ds)[2]],fDet,cutoff,outFile))
}

#=================================================================================================
# Steps 6 - 11 will be repeated 10 times for reporting each result and average
#                  (iSplitTimes = 10, default)
#=================================================================================================

source("s8.RegrrMethods.R")  # add external functions for regressions (all methods in one file!)

#-------------------------------------------------------------------------------------------------
# List with the results, with the same HEADER as the function are exporting
dfRes <- list("Regrrs" <- NULL, "Step" <- NULL, "tr.RMSE.10CV" <- NULL, "tr.R2.10CV" <- NULL,
              "tr.RMSESD.10CV" <- NULL, "tr.R2SD.10CV" <- NULL, "ts.RMSE.10CV" <- NULL, "ts.R2.10CV" <- NULL,
              "both.adjR2.10CV" <- NULL, "tr.adjR2.10CV" <- NULL, "ts.adjR2.10CV" <- NULL)

#-------------------------------------------------------------------------------------------------
for (i in 1:iSplitTimes) {                      # Step splitting number = i
  # -----------------------------------------------------------------------
  # (6) Dataset split: Training and Test sets
  # -----------------------------------------------------------------------
  cat("-> [6] Splitting dataset in Training and Test sets ...\n")
  
  source("s6.DsSplit.R")  # add external function
  iSeed=i                 # to reapeat the ds splitting, different values of seed will be used
  dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
  # get train and test from the resulted list
  ds.train<- dsList$train
  ds.test <- dsList$test
  
  # -----------------------------------------------------------------------
  # (7) Feature selection
  # -----------------------------------------------------------------------
  # Two types of regression funtions:
  # -> without wrapper methods (no feature selection excepts the built-in feature selection: Lasso, PLS)
  # -> with wrapper methods (all the functions without built-in feature selection)
  # ===>>>> if fFeatureSel = TRUE, the wrapper versions will be used in step 8
  
  # Note: additional feature selection could be implemented in the future
  
  # -----------------------------------------------------------------------
  # (8) Regressions
  # -----------------------------------------------------------------------
  #
  cat("-> [8] Run Regressions ...\n")
  
  # --------------------------------------------
  # 8.1. Basic LM : default
  # --------------------------------------------
  
  # -----------------------------------------------------------------------------------------
  # (8.2) GLM - based on AIC - Generalized Linear Model with Stepwise Feature Selection
  # -----------------------------------------------------------------------------------------
  if (fGLM==TRUE) {
    cat("-> [8.2] GLM stepwise - based on AIC ...\n")
    outFile <- file.path(PathDataSet,glmFile)   # the same folder as the input

    # both wrapper and nont-wrapper function in the same external file
    
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      my.stats<- GLMreg(ds.train,ds.test,fDet,outFile)   # run GLM
      my.stats$Step <- i                                 # add step number
      # print(data.frame(my.stats)) # print results for each split
      
    } else {                     # if there is a need for previous feature selection ->> use wrapper functions
      my.stats<- GLMregW(ds.train,ds.test,fDet,outFile)  # run GLM with wrapper method
    }
    
    #-------------------------------------------------------
    # Add output from GLM to the list of results
    #-------------------------------------------------------
    dfRes = mapply(c, my.stats, dfRes, SIMPLIFY=FALSE) # add results for one split to the main list of results
    
  }
  
  # --------------------------------------------
  # 8.3. PLS
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.4. Lasso
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.5. RBF
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.6. SVM radial
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.7. SVM linear
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.8. Neural Networks : default
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.9. SOM
  # --------------------------------------------
  
  # END OF REGRESSION Functions !!!
  
  #------------------------------------------------------------------------------
  # 9. Results from all models - ordered by 10 fold CV adjR2 (averaged values)
  #                             (if require, additional plots will be created)
  #-------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------
  # 10. Best model selection - detailed statistics
  #-------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------
  # 11. Test best model with test dataset
  #                   (+ Y randomization 100 times, bootstaping)
  #-------------------------------------------------------------------------------
  
}

#------------------------------------------------------------------------------
# 12. Report all results for 10 splittings
#-------------------------------------------------------------------------------
cat("[12] Results for all splitings\n")
print(data.frame(dfRes)) # print all results as data frame

#------------------------------------------------------------------------------
# 13. Assessment of Applicability Domain (plot leverage)
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 14. Output regression model as QMFR (to be implemented in the future)
#-------------------------------------------------------------------------------
