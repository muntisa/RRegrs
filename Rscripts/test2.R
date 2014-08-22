


#--------------------------------------------------------------------
# Write a LIST to CSV file
#--------------------------------------------------------------------
AppendList2CSv <- function(L,csvFile) {
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    write.table(names(l)[i], file=out_file, sep=",", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=FALSE)  #writes the name of the list elements ("A", "B", etc) 
    write.table(l[[i]], file=out_file, sep=",", dec=".", quote=FALSE, 
                col.names=NA, row.names=TRUE)  #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  


paramFile="Parameters.csv"
#-------------------------
# Read parameters from a file as data frame
# -----------------------------------------
Param.df <- read.csv(paramFile,header=T)

fDet         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fDet"),2])         # flag to calculate and print details for all the functions
fFilters     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fFilters"),2])     # flag to apply filters                          (2)
fScaling     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fScaling"),2])     # flag for dataset Scaling                       (3)
fRemNear0Var = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRemNear0Var"),2]) # flag for Removal of near zero variance columns (4)
fRemCorr     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRemCorr"),2])     # flag for Removal of correlated columns         (5)
fFeatureSel  = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fFeatureSel"),2])  # flag for wrapper methods for feature selection (7)

cutoff       = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="cutoff"),2]))  # cut off for correlated features
fLM          = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fLM"),2])     # flag to run LM            (8.1)
fGLM         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fGLM"),2])    # flag to run GLM           (8.2)
fPLS         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fPLS"),2])    # flag to run PLS           (8.3)
fLASSO       = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fLASSO"),2])  # flag to run LASSO         (8.4)
fRBFdda      = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRBFdda"),2]) # flat to run RBF DDA       (8.5)
fSVLM        = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fSVLM"),2])   # flat to run svmRadial     (8.6)
fNN          = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fNN"),2])     # flat to run NN            (8.8)

# ----------------------------------------------------------------------------------------
iScaling = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iScaling"),2])) # 1 = normalization; 2 = standardization, 3 = other; any other: no scaling
iScalCol = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iScalCol"),2])) # 1 = including dependent variable in scaling; 2: only all features; etc.
# ----------------------------------------------------------------------------------------
trainFrac   = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="trainFrac"),2]))   # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
iSplitTimes = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iSplitTimes"),2])) # default is 10; time to split the data in train and test (steps 6-11); report each step + average
noYrand     = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="noYrand"),2]))     # number of Y randomization (default = 100)

CVtypes = strsplit(as.character(Param.df[which(Param.df$RRegrs.Parameters=="CVtypes"),2]),";") # types of cross-validation methods

# -------------------------------------------------------------------------------------------------------
# Files
# -------------------------------------------------------------------------------------------------------
PathDataSet    = as.character(Param.df[which(Param.df$RRegrs.Parameters=="PathDataSet"),2]) # dataset folder for input and output files
DataFileName   = as.character(Param.df[which(Param.df$RRegrs.Parameters=="DataFileName"),2]) # input step 1 = ds original file name
No0NearVarFile = as.character(Param.df[which(Param.df$RRegrs.Parameters=="No0NearVarFile"),2]) # output step 3 = ds without zero near vars
ScaledFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ScaledFile"),2]) # output step 4 = scaled ds file name (in the same folder)
NoCorrFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="NoCorrFile"),2]) # output step 5 = dataset after correction removal

ResAvgs        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResAvgs"),2]) # the output file with averaged statistics for each regression method
ResBySplits    = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResBySplits"),2]) # the output file with statistics for each split and the averaged values
ResBest        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResBest"),2]) # the output file with statistics for the best model

lmFile         = as.character(Param.df[which(Param.df$RRegrs.Parameters=="lmFile"),2]) # LM output file for details
glmFile        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="glmFile"),2]) # GLM output file for details
plsFile        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="plsFile"),2]) # PLS output file for details
lassoFile      = as.character(Param.df[which(Param.df$RRegrs.Parameters=="lassoFile"),2]) # Lasoo Radial output file for details
rbfDDAFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="rbfDDAFile"),2]) # RBF DDA output file for details
svlmFile       = as.character(Param.df[which(Param.df$RRegrs.Parameters=="svlmFile"),2]) # SVM Radial output file for details
nnFile         = as.character(Param.df[which(Param.df$RRegrs.Parameters=="nnFile"),2]) # NN Radial output file for details

