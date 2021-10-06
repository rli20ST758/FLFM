

require(Rcpp)
require(RSpectra)
require(refund)
require(refund)
require(face)
require(splines)
require(MASS)
require(MFPCA)
require(mgcv)
require(Matrix)
require(fields)
require(rARPACK)
require(quadprog)
require(elasticnet)
require(pracma)
require(ggplot2)
require(gridExtra)
require(reshape2)

# Load functions
source("./R code/GenerateData.R")
source("./R code/GenerateData_L5.R")
source("./R code/FLFM.R")
source("./R code/Validation_L.R")
source("./R code/predict.FLFM.R")
source("./R code/predict.MFPCA.R")
source("./R code/fbps.cov.R") # for smooth covariance matrices


source("./R code/Validation_rho.R") # for sparse beta
source("./R code/spca_seqADMM.R") # for sparse beta
sourceCpp("./R code/eigendecomp.cpp") # for sparse beta
sourceCpp("./R code/MatrixMtp.cpp") # for sparse beta



## Generate data
Nsubj <- 200
Ntrain <- 100
Ntest <- 100
t <- 101
J <- 20
L_true <- 3 #L for generating data
FICC <- 2/3

set.seed(1)
data <- GenerateData(Nsubj=Nsubj, argval=t, J=J, L=L_true, FICC=FICC, SNR=5)
Y <- data$Y
Y_train <- lapply(1:J, function(x){(Y[[x]])[1:Ntrain,]})
Y_test <- lapply(1:J, function(x){(Y[[x]])[-(1:Ntrain),]})
# True values
beta_true <- as.matrix(data$beta)
C_true <- data$C
para <- data$lambda
G_true <- lapply(1:J, function(x){para[x] * data$G})
sigma2_true <- data$sigma2
argvals <- seq(0, 1, length=t)


#####################################################
## Get estimations by FLFM
#####################################################
# Set parameters
pve <- 0.95 #the proportion of variance explained 
knots <- 30 #the number of equidistant knots for all dimensions
p <- 3 #degrees of B-splines
m <- 2 #order of differencing penalty
tol <- 1e-4 #tolerance for iteration
Nmax <- 200  #number of maximal iteration

#####################################################
# Due to time limit, we put pound symbols on this part
# Select the best value L by cross validation. 
# seqL = c(2, 3, 4)
# L_result <- Validation_L(Y=Y_train, Nsubj=Ntrain, argvals=argvals, seqL=seqL, kfold=5, 
#                          knots=knots, p=p, m=m, tol=10e-4, Nmax=100)
# L <- L_result$L_min
# L_se <- L_result$L_min_se
# cv_diff <- L_result$cv_diff
#####################################################
L_se <- L_true

# Estimate parameters
time1 <- Sys.time()
FLFM_fit <- FLFM(Y=Y_train, Nsubj=Ntrain, argvals=argvals, L=L_se, pve1=pve, pve2=pve,
                 knots=knots, p=p, m=m, tol=tol, Nmax=Nmax, fit=T)
time2 <- Sys.time()
C_hat <-  FLFM_fit$C #estimated Cl
beta_hat <- FLFM_fit$beta #estimated beta
G_hat <- FLFM_fit$G #estimated Gj
sigma2_hat <- FLFM_fit$sigma2 # estimated sigma2

# Predict data in testing set
Pred <- predict.FLFM(object=FLFM_fit, newdata=Y_test, argvals=argvals)
Ypred_FLFM <- Pred$Ypred


#####################################################
## Get estimations by sparse FLFM
#####################################################
initial = list(C=FLFM_fit$C, beta=FLFM_fit$beta,
               G=FLFM_fit$G, sigma2=FLFM_fit$sigma2)

########################################################
# Due to time limit, we put pound symbols on this part
# Select the best value rho by cross validation
# seqrho <- seq(0, 0.5, 0.05) # FICC=2/3
# time2 <- Sys.time()
# rho_result <- Validation_rho(Y=Y_train, Nsubj=Ntrain, argvals=argvals, L=L_se, 
#                              seqrho=seqrho, kfold=5, initial=initial, knots=knots, 
#                              p=p, m=m, tol=10e-4, Nmax=100)
# rho <- rho_result$rho_min
########################################################
rho <- 0.2

# Estimate parameters
time1 <- Sys.time()
FLFM_sparse <- FLFM(Y=Y_train, Nsubj=Ntrain, argvals=argvals, L=L_se, sparse=TRUE, 
                    rho=rho, pve1=pve, pve2=pve, initial=initial, knots=knots, 
                    p=p, m=m, tol=tol, Nmax=Nmax,fit=TRUE)
time2 <- Sys.time()
sC_hat <-  FLFM_sparse$C #estimated Cl
sbeta_hat <- FLFM_sparse$beta #estimated beta
sG_hat <- FLFM_sparse$G #estimated Gj
ssigma2_hat <- FLFM_sparse$sigma2 # estimated sigma2

# Predict data in testing set
Pred_sparse <- predict.FLFM(FLFM_sparse, newdata=Y_test, argvals=argvals)
Ypred_sFLFM <- Pred_sparse$Ypred


#####################################################
## Get estimated results by MFPCA.
#####################################################
sim <- multiFunData(lapply(1:J, function(x){
  funData(argvals=argvals, X = Y_train[[x]])}))
# MFPCA based on univariate FPCA
M <- 80  
utype <- list(type="uFPCA", pve=pve)
time1 <- Sys.time()
MFPCA_fit <- MFPCA(sim, M=M, 
                   uniExpansions=list(utype,utype,utype,utype,utype,utype,utype,
                                      utype,utype,utype,utype,utype,utype,utype,
                                      utype,utype,utype,utype,utype,utype),fit=TRUE)
time2 <- Sys.time()
MFPCA_eigfun <- MFPCA_fit$functions
MFPCA_eigval <- MFPCA_fit$values

# Predict data in testing set
Pred_MFPCA <- predict.MFPCA(MFPCA_fit, data=Y_train, newdata=Y_test, 
                            argvals=argvals, pve=pve)
Ypred_MFPCA <- Pred_MFPCA$Ypred









