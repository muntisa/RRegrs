# -----------------------------------
# Example of RRegrs use
# -----------------------------------

# Libraries and external custom functions
# library(RRegrs)                 # load the RRegrs functions

# using the version in the package:
library(data.table)
library(corrplot)
source("../RRegrs/R/RRegrs_Functions.R")

# ------------------------------
# Default parameters of RRegrs
# ------------------------------
# DataFileName="ds.House.csv",PathDataSet="DataResults",
# ResAvgs="RRegsResAvgs.csv",ResBySplits="RRegrsResBySplit.csv",ResBest="RRegrsResBest.csv",

# fDet="T",fFilters="F",fScaling="T",fRemNear0Var="T",fRemCorr="T",fFeatureSel="F",

# fLM="T",fGLM="T",fPLS="T",fLASSO="T",fRBFdda="T",fSVRM="T",fNN="T",fRF="T",fSVMRFE="T",fENET="T",

# RFE_SVM_C="1;5;15;50",RFE_SVM_epsilon="0.01;0.1;0.3",
# cutoff="0.9",iScaling="1",iScalCol="1",trainFrac="0.75",iSplitTimes="2",noYrand="2",

# CVtypes="repeatedcv;LOOCV",

# No0NearVarFile="ds3.No0Var.csv",ScaledFile="ds4.scaled.csv",NoCorrFile="ds5.scaled.NoCorrs.csv",
# lmFile="8.1.LM.details.csv",glmFile="8.2.GLM.details.csv",plsFile="8.3.PLS.details.csv",
# lassoFile="8.4.LASSO.details.csv",rbfDDAFile="8.5.RBF_DDA.details.csv",svrmFile="8.6.SVMRadial.details.csv",
# nnFile="8.8.NN.details.csv",rfFile="8.9.RF.details.csv",svmrfeFile="8.10.SVMRFE.details.csv",
# enetFile="8.11.ENET.details.csv"

RRegrsResults = RRegrs(DataFileName="ds.House.csv",
                       fLM="T",fGLM="F",fPLS="T",fLASSO="F",fRBFdda="F",
                       fSVRM="F",fNN="F",fRF="T",fSVMRFE="F",fENET="F") # Run RRegrs without some methods
