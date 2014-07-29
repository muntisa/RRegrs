# ======================================================================
# Dataset spliting in Training and Test
# ======================================================================
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com
# ----------------------------------------------------------------------

# Inputs
# - ds = frame dataset object
# - fDet = flag for detais (TRUE/FALSE)
# - PathDataSet = pathway for results

# Output = training and test datasets
# if datails = TRUE, output files will be created

# ------------------------------------------
DsSplit <- function(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed) {
  my.datf<- ds
  # create TRAIN and TEST sets to build a model
  set.seed(iSeed)
  inTrain <- createDataPartition(1:dim(my.datf)[1],p = trainFrac,list = FALSE,groups=2)
							     # groups==2 forces to NOT partition
							     # based on quantiles of numeric values
  my.datf.train<- my.datf[inTrain,]            # TRAIN dataset frame         
  my.datf.test <- my.datf[-inTrain,]           # TEST dataset frame

  if (fDet == TRUE) {
    # write the TRAIN and TEST set files
    # the index of each row will in the dataset will not be saved (row.names=F)
    outTrain <- file.path(PathDataSet,"ds.Train.csv") # the same folder as the input
    write.csv(my.datf.train,outTrain,row.names=FALSE)
    outTest <- file.path(PathDataSet,"ds.Test.csv") # the same folder as the input
    write.csv(my.datf.test,outTest,row.names=FALSE) 
  }
  MyList<- list("train"=my.datf.train, "test"=my.datf.test) 
  return(MyList)  # return a list with training and test datasets
}
