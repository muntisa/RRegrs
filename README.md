RRegrs
======

eNanoMapper Developers |  [eNanoMapper Project] (http://www.enanomapper.net/)


The current tool is a collection of regression tools from R that could be used to search the best regression models for any dataset. The initial use of the script is aimed at finding QSAR models for chemoinformatics / nanotoxicology.

The full R script will contain: Loading dataset, Filter dataset, Scaling dataset, Feature selection, Regression models, Summary with top models, Statistics of the best model, etc.

The script will be modular in order to create flexible APIs.

The main authors are from the National Technical University of Athens (NTUA) and Maastricht University (UM) joining efforts for the EU eNanoMapper project.

The current implemented methods:
- Basic LM
- GLM based on AIC
- PLS
- Lasso
- RBF
- SVM radial
- Neural Networks

The methods to be implemented:
- SVM linear
- SOM
- Recursive Feature Extraction (SVM-RFE)

In addition, all methods will contain a wrapper version.

Outputs:
- CSV files for statistics
- PDF files for plots

For the best model, the last split of dataset will be used and additional files will be created.
The main statistics will be printed into:
- RRegrsResBySplit.csv = statistics for each split, regression method and cross-validation type
- RRegsResAvgs.csv     = averaged statistics for all splittings by method and cross-validation
- RRegrsResBest.csv    = best model statistics (best test adjR2)

To be done:
- Y randomization for the best model
- Applicability Domain for the best model 
- Corrections of formulas
- Wrapper functions
- Specific detailed outputs for each type of regression method
- General function for all regression methods with particular details for each regression type
