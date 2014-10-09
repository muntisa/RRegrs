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
# DataFileName="ds.House.csv",PathDataSet="DataResults",noCores=1,
# ResAvgs="RRegsResAvgs.csv",ResBySplits="RRegrsResAllSplits.csv",ResBest="RRegrsResBest.csv",
# fDet="T",fFilters="F",fScaling="T",fRemNear0Var="T",fRemCorr="T",fFeatureSel="F",
# fLM="T",fGLM="T",fPLS="T",fLASSO="T",fRBFdda="T",fSVRM="T",fNN="T",fRF="T",fSVMRFE="T",fENET="T",
# RFE_SVM_C="1;5;15;50",RFE_SVM_epsilon="0.01;0.1;0.3",
# cutoff=0.9,iScaling=1,iScalCol=1,trainFrac=0.75,iSplitTimes=10,noYrand=100,
# CVtypes="repeatedcv;LOOCV",
# No0NearVarFile="ds.No0Var.csv",ScaledFile="ds.scaled.csv",NoCorrFile="ds.scaled.NoCorrs.csv",
# lmFile="LM.details.csv",glmFile="GLM.details.csv",plsFile="PLS.details.csv",
# lassoFile="Lasso.details.csv",rbfDDAFile="RBF_DDA.details.csv",svrmFile="SVMRadial.details.csv",
# nnFile="NN.details.csv",rfFile="RF.details.csv",svmrfeFile="SVMRFE.details.csv",
# enetFile="ENET.details.csv"

# noCores = CPU cores: 0 = all CPU cores, 1 = no parallel, >1 = specific cores

RRegrsResults = RRegrs(DataFileName="ds.House.csv",noCores=2,iSplitTimes=2,noYrand=2,
                       fLM="T",fGLM="F",fPLS="T",fLASSO="T",fRBFdda="F",
                       fSVRM="F",fNN="F",fRF="T",fSVMRFE="F",fENET="F") # Run RRegrs without some methods
