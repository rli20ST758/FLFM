

predict.FLFM <- function(object, newdata, argvals){
  ############################################################ 
  #Arguments:
  #object -- an object returned by FLFM
  #newdata -- a data list consisting for which predicted values are desired. 
  #The list of length J contains densely observed multivariate (J-dimensional) functional data.
  #argvals -- the argument values of the function evaluations 
  ############################################################ 
  
  if(object$common_G == FALSE) {
    predict.LFM(object, newdata, argvals)
  } else {
    predict.LFM2(object, newdata, argvals)
  }
  
}



# FLFM using different individual covariance functions for functional outcomes
predict.LFM <- function(object, newdata, argvals){
  
  t <- length(argvals)
  Nsubj <- nrow(newdata[[1]])
  beta <- object$beta
  if(is.null(beta)) {
    L <- 0
    J <- length(object$G)
  } else{
    L <- ncol(beta)
    J <- nrow(beta)
  }
  lambda_C <- object$lambda_C
  lambda_G <- object$lambda_G
  phi <- object$phi
  psi <- object$psi
  sigma2 <- object$sigma2
  mu <- object$meanfunctions
  
  ####################################
  # Demean by fpca.face
  ####################################
  Y_resid <- list()
  for(j in 1:J){
    Y_resid[[j]] <- newdata[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
  }
  #Yi for each subject 
  Yi <- lapply(1:Nsubj, function(x){ 
    do.call(cbind, lapply(1:J, function(j){ Y_resid[[j]][x,] }) )
  }) 
  
  ####################################
  # Assign global eigenfunctions
  ####################################
  if(L == 0){
    lambda_G_temp = lambda_G
    psi_temp = psi
    # compute scores
    D <- lapply(1:J, function(x){
      if(length(lambda_G[[x]])==1) {
        lambda_G[[x]] * psi_temp[[x]]%*%t(psi_temp[[x]]) + diag(sigma2[x],t)
      } else{
        psi_temp[[x]]%*%diag(lambda_G[[x]])%*%t(psi_temp[[x]]) + diag(sigma2[x],t)
      } })
    invD <- lapply(1:J, function(x){solve(D[[x]])})
    zeta <- lapply(1:J, function(x){(lambda_G[[x]]*t(psi_temp[[x]]))%*%invD[[x]]%*%t(Y_resid[[x]])})
    Y_pred <- lapply(1:J, function(x){t(mu[x,] + psi_temp[[x]]%*%zeta[[x]])})
    xi = NULL
    zeta <- Reduce(rbind, zeta)
  }
  
  if(L > 0){
    # set NULL to 0
    lambda_C_temp = lambda_C
    lambda_G_temp = lambda_G
    phi_temp = phi
    psi_temp = psi
    for (x in 1:L) {
      if(length(lambda_C[[x]])==0){
        lambda_C_temp[[x]] = 1e-6
        phi_temp[[x]] = matrix(0, nc=1, nr=t)
      }
    }
    for (x in 1:J) {
      if(length(lambda_G[[x]])==0){
        lambda_G_temp[[x]] = 1e-6
        psi_temp[[x]] = matrix(0, nc=1, nr=t)
      }
    }
    
    # compute scores
    Mat1 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta[,x], phi_temp[[x]])}))
    Mat2 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta[,x]/sigma2, phi_temp[[x]])}))
    Mat3 <- lapply(1:J, function(x){t(psi_temp[[x]])%*%psi_temp[[x]]/sigma2[x]})
    Mat4 <- lapply(1:J, function(x){t(psi_temp[[x]])/sigma2[x]})
    invGamma1 <- bdiag(lapply(lambda_C_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
    invGamma2 <- bdiag(lapply(lambda_G_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
    A <- t(Mat1) %*% Mat2 + invGamma1
    C <- bdiag(Mat4) %*% Mat1
    D <- bdiag(Mat3) + invGamma2
    invD <- solve(D)
    MatE <- solve(A-t(C)%*%invD%*%C)
    MatF <- - MatE%*%t(C)%*%invD
    MatG <- - invD%*%C%*%MatE
    MatH <- invD - invD%*%C%*%MatF
    part_xi <- MatE%*%t(Mat2) + MatF%*%bdiag(Mat4)
    part_zeta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
    y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
    xi <- part_xi %*% y_vec
    zeta <- part_zeta %*% y_vec
    
    part1 <- Mat1%*%xi
    part2 <- bdiag(psi_temp) %*% zeta
    y_pred <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
    idx <- rep(1:J, Nsubj)
    Y_pred <- lapply(1:J, function(x){y_pred[idx==x,]})
    
    
    # delete scores corresponding to lambda=NULL
    idxC <- which(unlist(lambda_C_temp)==1e-6)
    if(length(idxC) > 0)  xi=as.matrix(xi[-idxC,])
    idxG <- which(unlist(lambda_G_temp)==1e-6)
    if(length(idxG) > 0)  zeta=as.matrix(zeta[-idxG,])
  }
  
  return(list(Ypred=Y_pred, xi=xi, zeta=zeta))
}



# FLFM using common individual covariance functions for all functional outcomes
predict.LFM2 <- function(object, newdata, argvals){
  
  t <- length(argvals)
  Nsubj <- nrow(newdata[[1]])
  beta <- object$beta
  L <- ncol(beta)
  J <- nrow(beta)
  lambda_C <- object$lambda_C
  lambda_G <- object$lambda_G
  phi <- object$phi
  psi <- object$psi
  gamma <- object$gamma
  sigma2 <- object$sigma2
  mu <- object$meanfunctions
  
  ####################################
  # Demean by fpca.face
  ####################################
  Y_resid <- list()
  for(j in 1:J){
    Y_resid[[j]] <- newdata[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
  }
  #Yi for each subject 
  Yi <- lapply(1:Nsubj, function(x){ 
    do.call(cbind, lapply(1:J, function(j){ Y_resid[[j]][x,] }) )
  }) 
  
  lambda_C_temp = lambda_C
  lambda_G_temp = lambda_G
  phi_temp = phi
  psi_temp = psi
  for (x in 1:L) {
    if(length(lambda_C[[x]])==0){
      lambda_C_temp[[x]] = 1e-6
      phi_temp[[x]] = matrix(0, nc=1, nr=t)
    }
  }
  
  # compute scores
  Mat1 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta[,x], phi_temp[[x]])}))
  Mat2 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta[,x]/sigma2, phi_temp[[x]])}))
  Mat3 <- lapply(1:J, function(x){gamma[x]^2*t(psi_temp)%*%psi_temp/sigma2[x]})
  Mat4 <- lapply(1:J, function(x){gamma[x]*t(psi_temp)/sigma2[x]})
  invGamma1 <- bdiag(lapply(lambda_C_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
  invGamma2 <- bdiag(lapply(1:J,function(x){if(length(lambda_G_temp)>1) diag(1/lambda_G_temp) else 1/lambda_G_temp}))
  A <- t(Mat1) %*% Mat2 + invGamma1
  C <- bdiag(Mat4) %*% Mat1
  D <- bdiag(Mat3) + invGamma2
  invD <- solve(D)
  MatE <- solve(A-t(C)%*%invD%*%C)
  MatF <- - MatE%*%t(C)%*%invD
  MatG <- - invD%*%C%*%MatE
  MatH <- invD - invD%*%C%*%MatF
  part_xi <- MatE%*%t(Mat2) + MatF%*%bdiag(Mat4)
  part_zeta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
  y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
  xi <- part_xi %*% y_vec
  zeta <- part_zeta %*% y_vec
  
  part1 <- Mat1%*%xi
  part2 <- bdiag(lapply(1:J, function(x){gamma[x]*psi_temp})) %*% zeta
  y_pred <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
  idx <- rep(1:J, Nsubj)
  Y_pred <- lapply(1:J, function(x){y_pred[idx==x,]})
  
  # delete scores corresponding to lambda=NULL
  idxC <- which(unlist(lambda_C_temp)==1e-6)
  if(length(idxC) > 0)  xi=as.matrix(xi[-idxC,])
  
  return(list(Ypred=Y_pred, xi=xi, zeta=zeta))
  
}



