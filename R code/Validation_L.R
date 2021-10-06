
Validation_L <- function(Y, Nsubj, argvals, seqL, sparse=FALSE, rho=0, kfold=5, 
                         knots=35, p=3, m=2, tol=10e-4, Nmax=100){
  ############################################################ 
  #Arguments:
  #Y -- observed data. A list of length J containing densely observed multivariate (J-dimensional) functional data. 
  #Nsubj -- number of subjects
  #argvals -- the argument values of the function evaluations 
  #seqL -- a sequence of L which are used in validation
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
  cv_diff <- matrix(0, nr=length(seqL), nc=kfold) #record the residual for each L
  est <- vector("list", length(seqL)*kfold)  #a list to save the estimated parameters
  
  for(L in seqL){
    print(paste("L is", L))
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
      LFM_result <- LFM_cv(Y=Y_train, Nsubj=Nsubj_train, argvals=argvals, L=L, sparse=sparse,
                           rho=rho, knots=knots, p=p, m=m, tol=tol, Nmax=Nmax)
      C_hat <-  LFM_result$C #estimated Cl
      beta_hat <- LFM_result$beta #estimated beta
      G_hat <- LFM_result$G #estimated Gj
      sigma2_hat <- LFM_result$sigma2 # estimated sigma2
      est[[(which(seqL==L)-1)*kfold+i_fold]] = LFM_result
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
          num_HS[j1,j2] <- sum(diff_Cov^2)/(t^2)
          den_HS[j1,j2] <- sum(Cov_sam^2)/(t^2)
        } #end j1
      } #end j2
      
     # Record the difference for each L and each fold
     #cv_diff[which(seqL==L), i_fold] <- sum(num_HS)/sum(den_HS) 
     #cv_diff[which(seqL==L), i_fold] <- sum(num_HS)
     cv_diff[which(seqL==L), i_fold] <- sum(num_HS)/(J^2)
     print(i_fold)
   } #end kfold
  } #end for every L
  
  ave <- rowSums(cv_diff)/kfold
  L_min <- seqL[which.min(ave)] 
  se <- sd(cv_diff[which.min(ave), ]) / sqrt(kfold)
  one_se <- min(ave) + se
  L_min_se <- seqL[ which(ave <= one_se)[1] ]
  initial <- est[[(which(seqL==L_min)-1)*kfold+which.min(cv_diff[which(seqL==L_min),])]]
  
  
  return( list(cv_diff=cv_diff, L_min = L_min, L_min_se = L_min_se, initial=initial))
  
}




LFM_cv <- function(Y, Nsubj, argvals, L, sparse=FALSE, rho=0, 
                   initial=list(C=NULL,beta=NULL,G=NULL,sigma2=NULL), 
                   knots=35, p=3, m=2, tol=1e-4, Nmax=100) {
  
  t <- length(argvals)
  J <- length(Y)
  ####################################
  ####step 1: Demean by fpca.face
  ####################################
  mu <- matrix(NA, nr=J, nc=t)
  sigma0 <- rep(NA, J)
  Y_resid <- list()
  for(j in 1:J){
    face_fit <- fpca.face(Y=Y[[j]], argvals=argvals, var=T)
    mu[j,] <- face_fit$mu
    # Demean for each variate
    Y_resid[[j]] <- Y[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
    sigma0[j] <- face_fit$sigma2
  }
  
  #Yi for each subject 
  Yi <- lapply(1:Nsubj, function(x){ 
    do.call(cbind, lapply(1:J, function(j){ Y_resid[[j]][x,] }) )
  }) 
  
  
  ############################################
  ####step 2: Estimate parameters by iteration
  ############################################
  # Initilize the values of parameters
  diff_C <- 1
  diff_beta <- 1
  diff_G <- 1
  diff_sigma2 <- 1
  if(!is.null(initial$C)) {
    C_new <- initial$C
  } else{
    C_new <- lapply(1:L, function(x){matrix(rnorm(t^2), nr=t, nc=t)})
  }
  if(!is.null(initial$beta)) {
    beta_new <- initial$beta
  } else{
    beta_new <- beta_init_func(Y=Y_resid, L=L)
  }
  if(!is.null(initial$G)) {
    G_new <- initial$G
  } else{
    G_new <- lapply(1:J, function(x){matrix(0, nr=t, nc=t)})
  }
  if(!is.null(initial$sigma2)) {
    sigma2_new <- initial$sigma2
  } else{
    sigma2_new <- sigma0
  }
  k = 0
  
  while( (diff_C>tol | diff_beta>tol | diff_G>tol | diff_sigma2>tol) & k <= Nmax){
    
    # Estimate each Cl
    C_old <- C_new
    C_new <- Renew_C(Yi=Yi, beta=beta_new, G=G_new, sigma2=sigma2_new, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    num_C <- unlist(lapply(1:L, function(x){norm(C_new[[x]]-C_old[[x]],"F")^2}))
    den_C <- unlist(lapply(1:L, function(x){norm(C_old[[x]],"F")^2}))
    diff_C <- sqrt(sum(num_C)/sum(den_C))
    
    # Estimate beta
    beta_old <- beta_new
    beta_new <- Renew_beta(Yi=Yi, C=C_new, beta=beta_old, G=G_new, sigma2=sigma2_new, Nsubj=Nsubj, L=L, J=J, sparse=sparse, rho=rho)
    diff_beta <- norm(beta_old-beta_new,"F")/norm(beta_old,"F")
    
    # Estimate each Gj and sigma2
    G_old <- G_new
    sigma2_old <- sigma2_new
    G_sigma <- Renew_G_sigma(Y_resid=Y_resid, C=C_new, beta=beta_new, sigma2=sigma2_old, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    G_new <- G_sigma$G
    sigma2_new <- G_sigma$sigma2
    num_G <- unlist(lapply(1:J, function(x){norm(G_new[[x]]-G_old[[x]],"F")^2}))
    den_G <- unlist(lapply(1:J, function(x){norm(G_old[[x]],"F")^2}))
    diff_G <- sqrt(sum(num_G)/sum(den_G))
    diff_sigma2 <- norm(sigma2_old-sigma2_new,"F")/norm(sigma2_old,"F")
    
    k = k + 1 
  } #end interation 
  
  
  ####################################
  ####step 3: Get estimated results 
  ####################################
  # Estimated Cl
  C_temp<- C_new 
  order <- order(unlist(lapply(1:L, function(x){sum(diag(C_temp[[x]]))})), decreasing=TRUE)
  C_hat <- lapply(1:L, function(x){C_temp[[order[x]]]})
  # Estimated beta
  beta_hat <- as.matrix(beta_new)
  beta_hat <- as.matrix(beta_hat[,order])
  # Estimated Gj
  G_hat <- G_new
  # Estimated sigma2
  sigma2_hat <- sigma2_new
  
  return(list(C=C_hat, beta=beta_hat, G=G_hat, sigma2=sigma2_hat))
}




