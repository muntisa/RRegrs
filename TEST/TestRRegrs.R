# -----------------------------------
# Example of RRegrs use
# -----------------------------------

# Libraries and external custom functions
# library(RRegrs)                 # load the RRegrs functions

# using the version in the package:
library(data.table)
library(corrplot)
source("../RRegrs/R/RRegrs_Functions.R")

paramFile = "Parameters.csv"  # set the file with parameters
RRegrsRes = RRegrs(paramFile) # Run RRegrs

