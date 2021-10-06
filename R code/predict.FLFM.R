
predict.FLFM <- function(object, newdata, argvals){
  ############################################################ 
  #Arguments:
  #object -- an object returned by LFM
  #newdata -- a data list consisting for which predicted values are desired. 
  #The list of length J contains densely observed multivariate (J-dimensional) functional data.
  #argvals -- the argument values of the function evaluations 
  ############################################################ 
  
  t <- length(argvals)
  Nsubj <- nrow(newdata[[1]])
  beta <- object$beta
  L <- ncol(beta)
  J <- nrow(beta)
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
  # efuncs1 <- list()
  # index <- 0
  # for(l in 1:L){
  #   efuncs1[[l]] <- kronecker(rep(1,J), phi[[l]])
  #   for(m in 1:ncol(efuncs1[[l]])){
  #     index <- index + 1
  #     assign(paste("Efunc1", index, sep=""), efuncs1[[l]][,m], envir=parent.frame()) #make it global variable
  #   }
  # }
  # 
  # max_n <- max(unlist(lapply(lambda_G, function(x)length(x))))
  # efuncs2 <- do.call(rbind, lapply(1:J, function(x){
  #   num <- length(lambda_G[[x]])
  #   if(num<max_n){return(cbind(psi[[x]],matrix(0,nr=t,nc=(max_n-num))))}
  #   return(psi[[x]])
  # }))
  # for(r in 1:max_n){
  #   assign(paste("Efunc2",r,sep=""), efuncs2[,r], envir = parent.frame()) #make it global variable
  # }
  # 
  # 
  # y_pred <- list()
  # xi <- matrix(NA, nrow=Nsubj, ncol=index)
  # eta0 <- matrix(NA, nrow=Nsubj, ncol=max_n*J)
  # g_var <- rep(as.factor(paste(1:J,sep="")),each=t) 
  # phimodel <- paste("s(g_var, by=Efunc1",1:index,", bs='re') ", collapse="+", sep="")
  # psimodel <- paste("s(g_var, by=Efunc2",1:max_n,", bs='re') ", collapse="+", sep="")
  # flags <- rep(1:L, unlist(lapply(1:L,function(x){length(lambda_C[[x]])})))
  # for (i in 1:Nsubj) {
  #   y <- as.vector(Yi[[i]])
  #   formula <- paste("y ~ ", phimodel, "+", psimodel)
  #   # Fit model and estimate scores
  #   fitting <- bam(as.formula(formula), drop.intercept=T)
  #   xi[i,] <- unlist(lapply(1:index, function(x){
  #     temp <- fitting$coefficients[((x-1)*J+1):(x*J)]
  #     sum(temp)/sum(beta_hat[,flags[x]])
  #   }))
  #   eta0[i,] <- fitting$coefficients[-(1:(index*J))]
  #   y_pred[[i]] <- matrix(fitting$fitted.values, nc=t, byrow=T) + mu
  # }
  # 
  # Y_pred <- lapply(1:J, function(x){do.call("rbind",lapply(1:Nsubj, function(sub){y_pred[[sub]][x,]}))})
  # eta <- lapply(1:J, function(x){
  #   start <- max_n*(x-1)+1
  #   end <- max_n*(x-1) + length(lambda_G[[x]])
  #   eta0[,start:end]})
  
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
  part_eta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
  y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
  xi <- part_xi %*% y_vec
  eta <- part_eta %*% y_vec
  
  part1 <- Mat1%*%xi
  part2 <- bdiag(psi_temp) %*% eta
  y_pred <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
  idx <- rep(1:J, Nsubj)
  Y_pred <- lapply(1:J, function(x){y_pred[idx==x,]})


  # delete scores corresponding to lambda=NULL
  idxC <- which(unlist(lambda_C_temp)==1e-6)
  if(length(idxC) > 0)  xi=as.matrix(xi[-idxC,])
  idxG <- which(unlist(lambda_G_temp)==1e-6)
  if(length(idxG) > 0)  eta=as.matrix(eta[-idxG,])
  
  
  return(list(Ypred=Y_pred, xi=xi, eta=eta))
  
}


