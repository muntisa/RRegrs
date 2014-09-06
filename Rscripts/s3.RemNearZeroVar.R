#=======================================================================
# Removal of near zero variance columns
#=======================================================================
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com
#-----------------------------------------------------------------------
# inputs:
# - ds = dataset frame
# - fDet = flag for detais (TRUE/FALSE)
# - outFileName = new file name  (it could include the path)

# output = ds.Rem0NearVar  (ds without columns with near zero variance)
# if datails = TRUE, output the new ds as a file
# ------------------------------------------

RemNear0VarCols <- function(ds,fDet=FALSE,outFile="ds3.No0Var.csv") {
  # default parameters are no details, with a CSV file name
  library(caret)
  ds.Rem0NearVar <- ds               # default output without any modification
  ds.var <- nearZeroVar(ds)          # get the near zero columns
  if (!length(ds.var) == FALSE) {    # remove the columns only if nearZeroVar identified; if no columns to remove, ds will be the same
    ds.Rem0NearVar <- ds[,-(ds.var)] # get only the columns without this problem
    if (fDet == TRUE) {              # write as details the corrected ds file
      write.csv(ds.Rem0NearVar, outFile,row.names=F, quote=F)
    }
  }
  return(as.data.frame(ds.Rem0NearVar)) # return the new data frame without near zero variance
}
