
Validation_rho <- function(Y, Nsubj, argvals, L, seqrho, kfold=5, 
                           initial=list(C=NULL,beta=NULL,G=NULL,sigma2=NULL),
                           knots=35, p=3, m=2, tol=10e-4, Nmax=100){
  ############################################################ 
  #Arguments:
  #Y -- observed data. A list of length J containing densely observed multivariate (J-dimensional) functional data. 
  #Nsubj -- number of subjects
  #argvals -- the argument values of the function evaluations 
  #seqrho -- a sequence of rho which are used in validation
  #initial --  the initial values of parameters for the iteration. It can be from cross-validation for L. 
  #            It highly recommends to set initial values as the estimation for sparse=FALSE.
  #knots -- list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #tol -- tolerance for iteration, default to 10e-4
  #Nmax -- the number of maximal iteration, default to 100
  ############################################################
  
  #############################
  # Demean by fpca.face
  #############################
  t <- length(argvals)
  J <- length(Y)
  mu <- matrix(NA, nr=J, nc=t)
  Y_resid <- list()
  for(j in 1:J){
    face_fit <- fpca.face(Y=Y[[j]], argvals=argvals, var=T)
    mu[j,] <- face_fit$mu
    # Demean for each variate
    Y_resid[[j]] <- Y[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
  }
  
  ID <- 1:Nsubj
  folds <- cv.folds(Nsubj, folds=kfold) #index for k fold
  cv_diff <- matrix(0, nr=length(seqrho), nc=kfold) #record the residual for each L
  
  for(rho in seqrho){
    print(paste("The rho is", rho))
    for(i_fold in 1:kfold){
      
      #Set the training set
      ID_train <- which(!(ID %in% folds[[i_fold]]))
      Nsubj_train <- length(ID_train)
      Y_train <- lapply(1:J, function(x){Y_resid[[x]][ID_train, ] })
      #Set the validation set
      ID_test <- ID[-ID_train]
      Nsubj_test <- length(ID_test)
      Y_test <- lapply(1:J, function(x){Y_resid[[x]][ID_test, ]})
      
      #############################
      # Estimate parameters
      #############################
      LFM_result <- LFM_cv(Y=Y_train, Nsubj=Nsubj_train, argvals=argvals, L=L, sparse=TRUE, 
                           rho=rho, initial=initial, knots=knots, p=p, m=m, tol=tol, Nmax=Nmax)
      C_hat <-  LFM_result$C #estimated Cl
      beta_hat <- LFM_result$beta #estimated beta
      G_hat <- LFM_result$G #estimated Gj
      sigma2_hat <- LFM_result$sigma2 # estimated sigma2
      
      ###################################################################################
      # Compute the difference between validation estimate and the predict covariance
      ###################################################################################
      num_HS <- matrix(0, nrow=J, ncol=J)
      den_HS <- matrix(0, nrow=J, ncol=J)
      for(j1 in 1:J){
        for(j2 in 1:J){
          # Get sample corvariance for data in validation set
          if(j1==j2) {
            Cov_sam <- t(Y_test[[j1]]) %*% Y_test[[j2]]/Nsubj_test - diag(sigma2_hat[j1],t)
            Cov_pre <- Reduce('+', lapply(1:L, function(x){beta_hat[j1,x]*beta_hat[j2,x]*C_hat[[x]]}))
            Cov_pre  <- Cov_pre  + G_hat[[j1]]
          }else{
            Cov_sam <- t(Y_test[[j1]]) %*% Y_test[[j2]] / Nsubj_test
            Cov_pre <- Reduce('+', lapply(1:L, function(x){beta_hat[j1,x]*beta_hat[j2,x]*C_hat[[x]]}))
          }
          
          diff_Cov <- Cov_sam  - Cov_pre
          num_HS[j1,j2] <- sum(diff_Cov^2)
          den_HS[j1,j2] <- sum(Cov_sam^2)
        } #end j1
      } #end j2
      
      # Record the difference for each L and each fold
      #cv_diff[which(seqrho==rho), i_fold] <- sum(num_HS)/sum(den_HS) 
      cv_diff[which(seqrho==rho), i_fold] <- sum(num_HS)/(J^2*t^2)
      print(i_fold)
    } #end kfold
  } #end for every rho
  
  ave <- rowSums(cv_diff)/kfold
  rho_min <- seqrho[which.min(ave)] 
  se <- sd(cv_diff[which.min(ave), ]) / sqrt(kfold)
  one_se <- min(ave) + se
  rho_min_se <- seqrho[ which(ave <= one_se)[1] ]
  
  return( list(cv_diff=cv_diff, rho_min = rho_min, rho_min_se = rho_min_se))
  
}





