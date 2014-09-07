#===================================================================================
# Scaling dataset (Step 4)
#===================================================================================
# s = { 1,2,3 } - type of scaling: 1 = normalization, 2 = standardization, 3 = other
# c = the number of column into the dataset to start scaling
# - if c = 1: included the dependent variable
# - if c = 2: only the features will be scaled

# fDet = if details need to be printed    (TRUE/FALSE)
# inFileName  = file name      (it could include the path)
# outFileName = new file name  (it could include the path)
#-----------------------------------------------------------------------------------

ScalingDS <- function(ds,s=1,c=1,fDet=FALSE,outFileName="ds4.scaled.csv") {
  # Default scaling = NORMALIZATION !
  # if we use preProcess() from caret, errors are generating if threre are zero variance columns!
  
  # DEFAULT scaled dataset = original
  # if other s diffent of 1,2,3 is used, no scaling!
  DataSet.scaled <- ds
  
  # if NORMALIZATION
  if (s==1) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- ((ds-min(ds))/(max(ds)-min(ds))) # normalize all the columns
  }
  # if STADARDIZATION
  if (s==2) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- scale(ds[c:ncol(ds)],center=TRUE,scale=TRUE)  
  }
  
  # if other scaling
  if (s==3) {
    # Scale all the features (from feature 2 bacause feature 1 is the predictor output)
    # TO ADD THE CODE ! 
  }
  
  # if DETAILS
  if (fDet ==TRUE) {
    # write the result into a separated file
    write.csv(DataSet.scaled, outFileName,row.names=F, quote=F)  
  }
  return (as.data.frame(DataSet.scaled)) # return the scaled data frame
}
