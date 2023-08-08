# FLFM
This project provides supplementary materials for the paper, Latent Factor Model for Multivariate Functional Data. 
- The proposed method is implemented in the function FLFM() in the folder "R code".
- An example to run FLFM() for a simulated data is shown in the file demo.R.
- An example to run FLFM() for an EEG data is shown in the file demo_EEG.R. The EEG data is saved in the folder data.
- The RMD file provides examples of how to use the function FLFM(). Moreover, comparisons with MFPCA proposed in Happ and Greven (2018) are evaluated. It would help understand the proposed method. 
To run the R function FLFM(), R packages "refund", "face", "splines", "MASS", "mgcv", "Matrix", "fields", "rARPACK", "glmnet", "quadprog", "elasticnet", "pracma", "Rcpp" and "RSpectra" are required.