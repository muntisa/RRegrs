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

library(caret)

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
fLM          = TRUE  # flag to run LM            (8.1)
fGLM         = TRUE  # flag to run GLM           (8.2)
fPLS         = FALSE  # flag to run PLS           (8.3)
fLASSO       = FALSE  # flag to run LASSO         (8.4)
fRBFdda      = FALSE  # flat to run RBF DDA       (8.5)
fSVLM        = FALSE # flat to run svmRadial.RMSE (8.6)
fNN          = FALSE  # flat to run NN            (8.8)

# ----------------------------------------------------------------------------------------
iScaling = 1 # 1 = normalization; 2 = standardization, 3 = other; any other: no scaling
iScalCol = 1 # 1 = including dependent variable in scaling; 2: only all features; etc.
# ----------------------------------------------------------------------------------------
trainFrac   = 3/4 # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
iSplitTimes = 2   # default is 10; time to split the data in train and test (steps 6-11); report each step + average
noYrand     = 3   # number of Y randomization (default = 100)

CVtypes    <- c("repeatedcv","LOOCV")             # cross-validation types: 10-CV and LOOCV

# -------------------------------------------------------------------------------------------------------
# Files
# -------------------------------------------------------------------------------------------------------
PathDataSet    = "DataResults"            # dataset folder for input and output files
DataFileName   = "ds.csv"                 # input step 1 = ds original file name
No0NearVarFile = "ds3.No0Var.csv"         # output step 3 = ds without zero near vars
ScaledFile     = "ds4.scaled.csv"         # output step 4 = scaled ds file name (in the same folder)
NoCorrFile     = "ds5.scaled.NoCorrs.csv" # output step 5 = dataset after correction removal

ResAvgs        = "RRegsResAvgs.csv"       # the output file with averaged statistics for each regression method
ResBySplits    = "RRegrsResBySplit.csv"   # the output file with statistics for each split and the averaged values
ResBest        = "RRegrsResBest.csv"      # the output file with statistics for the best model

lmFile         = "8.1.LM.details.csv"           # LM output file for details
glmFile        = "8.2.GLM.details.csv"          # GLM output file for details
plsFile        = "8.3.PLS.details.csv"          # PLS output file for details
lassoFile      = "8.4.LASSO.details.csv"        # Lasoo Radial output file for details
rbfDDAFile     = "8.5.RBF_DDA.details.csv"      # RBF DDA output file for details
svlmFile       = "8.6.SVMRadial.details.csv"    # SVM Radial output file for details
nnFile         = "8.8.NN.details.csv"           # NN Radial output file for details

# Generate path + file name = original dataset
inFile <- file.path(PathDataSet, DataFileName)

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
if (fFilters==TRUE) {
  # cat("-> [2] Filtering dataset ... \n")
}

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
# Initialize the list with the statistics results; the same HEADER as the function output
dfRes <- list("RegrMeth"     = NULL,
              "Split No"     = NULL,    
              "CVtype"       = NULL,      
              "NoModelFeats" = NULL,
              "ModelFeats"   = NULL,
              "adjR2.tr"  = NULL,
              "RMSE.tr"   = NULL,
              "R2.tr"     = NULL,
              "RMSEsd.tr" = NULL,
              "R2sd.tr"   = NULL,
              "adjR2.ts"= NULL,
              "RMSE.ts" = NULL,
              "R2.ts"   = NULL,
              "corP.ts" = NULL,
              "adjR2.both" = NULL,
              "RMSE.both"  = NULL,
              "R2.both"    = NULL)

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
  
  # Training data set info - to be added later!
  # ------------------------------------------------------------------------
  #nCases     <- dim(ds)[1]                     # number of total cases
  #nCases.tr  <- dim(ds.train)[1]               # number of training cases
  #nFeatures  <- dim(ds.train)[2] - 1           # number of input features (it could be different with the fitted model features!)
  #FeatureList<- paste(names(ds.train)[2:nFeatures], collapse="+") # list of input feature names, from second name, first name is the output variable
  #OutVar     <- names(ds.train)[1]             # name of predicted variable (dependent variable)
  
  # -----------------------------------------------------------------------
  # (7) Feature selection
  # -----------------------------------------------------------------------
  # Two types of regression funtions:
  # -> without wrapper methods (no feature selection excepts the built-in feature selection: Lasso, PLS)
  # -> with wrapper methods (all the functions without built-in feature selection)
  # ===>>>> if fFeatureSel = TRUE, the wrapper versions will be used in step 8
  
  # Note: additional feature selection could be implemented in the future
  
  # -----------------------------------------------------------------------
  # (8) REGRESSION METHODS
  # -----------------------------------------------------------------------
  #
  cat("-> [8] Run Regressions ...\n")
  
  # --------------------------------------------
  # 8.1. Basic LM : default
  # --------------------------------------------
  if (fLM==TRUE) {   # if LM was selected, run the method
    cat("-> [8.1] LM ...\n")
    outFile.LM <- file.path(PathDataSet,lmFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) { # there is no CV but it will be implemented in the future!!!
        my.stats.LM   <- LMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LM) # run GLM for each CV and regr method
        
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.LM, dfRes, SIMPLIFY=FALSE)
        
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      # run LM with wrapper method (TO BE IMPLEMENTED!)
    }
    
  } # end LM
  
  # -----------------------------------------------------------------------------------------
  # (8.2) GLM based on AIC regression - Generalized Linear Model with Stepwise Feature Selection
  # -----------------------------------------------------------------------------------------
  if (fGLM==TRUE) {   # if GLM was selected, run the method
    cat("-> [8.2] GLM stepwise - based on AIC ...\n")
    outFile.GLM <- file.path(PathDataSet,glmFile)   # the same folder as the input is used for the output

    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      
      # List with data set and method information
      #my.stats.dsInfo <- list("RegrMethod"= RegrMethod,           # regression method name
                              #"NoCases"= as.numeric(nCases),      # no. of cases
                              #"InNoVars"= as.numeric(nFeatures),  # no. of input features
                              #"InFeatures"= FeatureList,          # list with feature names
                              #"PredVar" = OutVar)                 # name of predicted variable (dependent variable)
      
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.GLM   <- GLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.GLM) # run GLM for each CV and regr method
        #my.stats.split <- c(my.stats.dsInfo,my.stats.GLM) # merge the ds info with statistics results for each Cv & reg method
        
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.GLM, dfRes, SIMPLIFY=FALSE)
      
        } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
      # run GLM with wrapper method (TO BE IMPLEMENTED!)
    }
  
  } # end GLM
  
  # --------------------------------------------
  # 8.3. PLS
  # --------------------------------------------
  if (fPLS==TRUE) {   # if LASSO was selected, run the method
    outFile.PLS <- file.path(PathDataSet,plsFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      cat("-> [8.3] PLS ...\n")
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.PLS  <- PLSreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run SVLM Radial for each CV and regr method
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=FALSE)
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      cat("-> [8.3] PLS Wrapper Feature Selection ...\n")
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.PLS  <- PLSregWSel(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run SVLM Radial for each CV and regr method
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=FALSE)
      } # end CV types 
    }
    
  } # end PLS with wrapper
  
  # --------------------------------------------
  # 8.4. LASSO regression
  # --------------------------------------------
  if (fLASSO==TRUE) {   # if LASSO was selected, run the method
    cat("-> [8.4] LASSO ...\n")
    outFile.LASSO <- file.path(PathDataSet,lassoFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.LASSO  <- LASSOreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LASSO) # run SVLM Radial for each CV and regr method
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.LASSO, dfRes, SIMPLIFY=FALSE)
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      # run Lasso with wrapper method (TO BE IMPLEMENTED!)
    }
    
  } # end Lasso
  
  # --------------------------------------------
  # 8.5. RBF network with the DDA algorithm regression (caret)
  # --------------------------------------------
  if (fRBFdda==TRUE) {   # if SVM Radial was selected, run the method
    cat("-> [8.6] RBF network with the DDA ...\n")
    outFile.rbfDDA <- file.path(PathDataSet,rbfDDAFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.rbfDDA  <- RBF_DDAreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.rbfDDA) # run SVLM Radial for each CV and regr method
        #-------------------------------------------------------
        # Add output from SVM Radial to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.rbfDDA, dfRes, SIMPLIFY=FALSE)
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      # run rbfDDA with wrapper method (TO BE IMPLEMENTED!)
    }
    
  } # end rbfDDA
  
  # --------------------------------------------
  # 8.6. SVM radial regression
  # --------------------------------------------
  if (fSVLM==TRUE) {   # if SVM Radial was selected, run the method
    cat("-> [8.6] SVM radial ...\n")
    outFile.SVLM <- file.path(PathDataSet,svlmFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.SVLM  <- SVLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVLM) # run SVLM Radial for each CV and regr method
        #-------------------------------------------------------
        # Add output from SVM Radial to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.SVLM, dfRes, SIMPLIFY=FALSE)
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      # run SVLM with wrapper method (TO BE IMPLEMENTED!)
    }
    
  } # end SVLM
  
  # --------------------------------------------
  # 8.7. SVM linear
  # --------------------------------------------
  
  # --------------------------------------------
  # 8.8. Neural Networks Regression
  # --------------------------------------------
  if (fNN==TRUE) {   # if NNet was selected, run the method
    cat("-> [8.8] Neural Networks ...\n")
    outFile.NN <- file.path(PathDataSet,nnFile)   # the same folder as the input is used for the output
    
    # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
    if (fFeatureSel==FALSE) {    # if there is no need of feature selection ->> use normal functions
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) {
        my.stats.NN  <- NNreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.NN) # run NNet for each CV and regr method
        #-------------------------------------------------------
        # Add output from NNet to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.NN, dfRes, SIMPLIFY=FALSE)
      } # end CV types
    } 
    else    # if there is a need for previous feature selection ->> use wrapper functions
    {                     
      # run NNet with wrapper method (TO BE IMPLEMENTED!)
    }
    
  } # end NNet
  # --------------------------------------------
  # 8.9. SOM
  # --------------------------------------------
  
  # END OF REGRESSION Functions !!!
}

#------------------------------------------------------------------------------
# 9. Results for all splittings (not ordered)
#-------------------------------------------------------------------------------
cat("[12] Results for all splitings ...\n")
df.res <- data.frame(dfRes)
print(df.res) # print all results as data frame

# Writing the statistics into output files: one with detailed splits, other with only averages
# File names includin paths for the statistics outputs (only averages and split detailed; +averages [to be implemented])
ResBySplitsF <- file.path(PathDataSet,ResBySplits)   # the main output file with statistics for each split
write.csv(df.res, file = ResBySplitsF)    # write statistics data frame into a CSV output file
# file.show(ResBySplitsF)  # show the statistics file!

#-------------------------------------------------------------------------------------
# Averaged values of the results by each Regression Method & CV type
#-------------------------------------------------------------------------------------
cat("-> Averaged statistics ...\n")

ResAvgsF <- file.path(PathDataSet,ResAvgs)           # the main output file with averaged statistics for each regression method
library(data.table)
dt.res  <- data.table(df.res) # convert data frame into data table (for sorting abd averaging)

# MEANS for each Regression Method & CV type
#--------------------------------------------------------------------------------------------------------------
# means for all CV types, not only 10CV
dt.mean <- dt.res[,list(adjR2.tr.Avg=mean(adjR2.tr),RMSE.tr.Avg=mean(RMSE.tr),R2.tr.Avg=mean(R2.tr),
                        RMSEsd.tr.Avg=mean(RMSEsd.tr),R2sd.tr.Avg=mean(R2sd.tr),adjR2.ts.Avg=mean(adjR2.ts),
                        RMSE.ts.Avg=mean(RMSE.ts),R2.ts.Avg=mean(R2.ts),corP.ts.Avg=mean(corP.ts),
                        adjR2.both.Avg=mean(adjR2.both),RMSE.both.Avg=mean(RMSE.both),
                        R2.both.Avg=mean(R2.both)),by="RegrMeth,CVtype"]

dt.mean     <- dt.mean[dt.mean$CVtype=="repeatedcv",]   # keep only the 10CV results to be used to find the best model
dt.mean.ord <- dt.mean[order(-rank(adjR2.ts.Avg))]      # descendent order the averages by adjR2.ts.Avg

# Write averages descendent ordered by adjR2.ts.Avg
#-------------------------------------------------------------------------------
write.csv(data.frame(dt.mean.ord), file = ResAvgsF)    # write statistics data frame into a CSV output file

#------------------------------------------------------------------------------
# 10. Best model selection - detailed statistics
#-------------------------------------------------------------------------------
cat("-> Best model analysis ...\n")

# ADD an algorithm to verifty similar adjR2 values:
# From the best ones (+/- 0.05 of adjR2), chose the one with less variables, after that the one with min RMSE!!!

best.dt  <- dt.mean.ord[1] # the best model should be the first value in the descendent ordered results
best.reg <- paste(best.dt$RegrMeth,collapse="") # best regrression method

#----------------------------------------------------------
# 11. Best model detailed statistics 
#----------------------------------------------------------
# Write the best model statistics
ResBestF <- file.path(PathDataSet,ResBest)
write.table("Averaged values for all spits: ", file = ResBestF, append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
# write.csv(data.frame(best.dt), file = ResBestF)    # write statistics data frame into a CSV output file
write.table(data.frame(best.dt), file = ResBestF,append = TRUE, sep = ",",col.names = TRUE,quote = FALSE) # write statistics data frame into a CSV output file

# Use the last split for dataset (ds.train & ds.test) ! (or chose other one?)

# Run the caret function with the method from the best method, for one training-test split only
# and append the details in the best model output file

if (best.reg=="lm") {
  my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run GLM for each CV and regr method
}
if (best.reg=="glmStepAIC") {
  my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run GLM for each CV and regr method
}
if (best.reg=="pls") {
  my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run SVLM Radial for each CV and regr method
}
if (best.reg=="lasso") {
  my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run SVLM Radial for each CV and regr method
}
if (best.reg=="rbfDDA") {  
  my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run SVLM Radial for each CV and regr method
}
if (best.reg=="svmRadial") {  
  my.stats.reg  <- SVLMreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run SVLM Radial for each CV and regr method
}
if (best.reg=="nnet") {  
  my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,TRUE,ResBestF) # run NNet for each CV and regr method
} 

#------------------------------------------------------------------------------
# 12. Test best model with test dataset
#                   (+ Y randomization 100 times, bootstaping)
#-------------------------------------------------------------------------------
# ratios Yrand R2 - Best model R2 / Best model R2
R2Diff.Yrand <- Yrandom(best.reg,my.stats.reg$R2.ts,noYrand,ResBestF) 

#------------------------------------------------------------------------------
# 13. Assessment of Applicability Domain (plot leverage)
#-------------------------------------------------------------------------------

#----------------------------
# Indicate main result files
#----------------------------
cat("\n MAIN RESULT FILES:\n")
cat("===================================\n")
cat("Statistics for each data set splitting/method/CV type:",ResBySplits,"\n")
cat("Averages for all data set splittings by method/CV type:",ResAvgsF,"\n")
cat("Best model statistics:",ResBestF,"\n")
cat("\n* if you choose Details, additional CSV files will be create for each method.\n")

